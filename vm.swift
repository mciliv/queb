import Cocoa
import Virtualization

// MARK: - Helpers

func expandedURL(from path: String) -> URL {
    URL(fileURLWithPath: (path as NSString).expandingTildeInPath)
}

let supportDir = URL(fileURLWithPath: NSHomeDirectory())
    .appendingPathComponent("Library/Application Support/UbuntuVM", isDirectory: true)

// MARK: - CLI

let args = CommandLine.arguments
let command = args.count > 1 ? args[1] : "start"

func printHelp() {
    print("""
    ubuntu-vm — Ubuntu ARM64 VM using Apple Virtualization.framework

    USAGE
      ubuntu-vm [command]

    COMMANDS
      start     Launch the VM in a GUI window (default)
      restart   Kill running VM, wipe disk/EFI state, start fresh
      -h        Show this help

    DETAILS
      CPU/RAM   4 cores, 4 GB
      Disk      20 GB sparse image at ~/Library/Application Support/UbuntuVM/disk.img
      EFI       ~/Library/Application Support/UbuntuVM/efi-variable-store
      ISO       Auto-detects ~/Downloads/ubuntu-25.10-desktop-arm64.iso
                or first .iso in ~/Downloads
      Network   NAT (appears as wired ethernet in Ubuntu, not WiFi)

    BUILD
      swiftc vm.swift -o ubuntu-vm -framework Virtualization -framework Cocoa
      codesign --force --sign - --entitlements entitlements.plist ubuntu-vm
    """)
}

switch command {
case "-h", "--help", "help":
    printHelp()
    exit(0)

case "restart":
    let kill = Process()
    kill.executableURL = URL(fileURLWithPath: "/usr/bin/pkill")
    kill.arguments = ["-f", "UbuntuVM"]
    try? kill.run()
    kill.waitUntilExit()
    try? FileManager.default.removeItem(at: supportDir)
    print("Cleared VM state. Starting fresh...")
    fallthrough

case "start":
    break

default:
    printHelp()
    exit(0)
}

// MARK: - VM Configuration

func createVMConfiguration() throws -> VZVirtualMachineConfiguration {
    let config = VZVirtualMachineConfiguration()
    config.cpuCount = 4
    config.memorySize = 4 * 1024 * 1024 * 1024

    // Audio
    let outputStream = VZVirtioSoundDeviceOutputStreamConfiguration()
    outputStream.sink = VZHostAudioOutputStreamSink()
    let outputDevice = VZVirtioSoundDeviceConfiguration()
    outputDevice.streams = [outputStream]

    let inputStream = VZVirtioSoundDeviceInputStreamConfiguration()
    inputStream.source = VZHostAudioInputStreamSource()
    let inputDevice = VZVirtioSoundDeviceConfiguration()
    inputDevice.streams = [inputStream]

    config.audioDevices = [outputDevice, inputDevice]

    // ISO
    let preferredISOPath = "~/Downloads/ubuntu-25.10-desktop-arm64.iso"
    var ubuntuISOURL = expandedURL(from: preferredISOPath)

    if !FileManager.default.fileExists(atPath: ubuntuISOURL.path) {
        let downloadsURL = expandedURL(from: "~/Downloads")
        if let iso = try? FileManager.default.contentsOfDirectory(at: downloadsURL, includingPropertiesForKeys: nil)
            .first(where: { $0.pathExtension.lowercased() == "iso" }) {
            ubuntuISOURL = iso
        }
    }

    guard FileManager.default.fileExists(atPath: ubuntuISOURL.path) else {
        fatalError("Ubuntu ISO not found. Download arm64 desktop ISO to ~/Downloads.")
    }

    // Disk image
    try? FileManager.default.createDirectory(at: supportDir, withIntermediateDirectories: true)
    let diskImageURL = supportDir.appendingPathComponent("disk.img")
    if !FileManager.default.fileExists(atPath: diskImageURL.path) {
        let fd = open(diskImageURL.path, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR)
        guard fd != -1 else { fatalError("Failed to create disk image") }
        let size: off_t = 20 * 1024 * 1024 * 1024
        guard ftruncate(fd, size) == 0 else { close(fd); fatalError("Failed to size disk") }
        close(fd)
    }

    let writableAttachment = try VZDiskImageStorageDeviceAttachment(url: diskImageURL, readOnly: false)
    let isoAttachment = try VZDiskImageStorageDeviceAttachment(url: ubuntuISOURL, readOnly: true)
    config.storageDevices = [
        VZVirtioBlockDeviceConfiguration(attachment: writableAttachment),
        VZVirtioBlockDeviceConfiguration(attachment: isoAttachment)
    ]

    // EFI boot loader
    let bootLoader = VZEFIBootLoader()
    let efiStore = supportDir.appendingPathComponent("efi-variable-store")
    if !FileManager.default.fileExists(atPath: efiStore.path) {
        bootLoader.variableStore = try VZEFIVariableStore(creatingVariableStoreAt: efiStore)
    } else {
        bootLoader.variableStore = VZEFIVariableStore(url: efiStore)
    }
    config.bootLoader = bootLoader
    config.platform = VZGenericPlatformConfiguration()

    // Graphics
    let graphicsDevice = VZVirtioGraphicsDeviceConfiguration()
    graphicsDevice.scanouts = [VZVirtioGraphicsScanoutConfiguration(widthInPixels: 1920, heightInPixels: 1080)]
    config.graphicsDevices = [graphicsDevice]

    // Network (NAT — shows as wired ethernet in Ubuntu, not WiFi)
    let net = VZVirtioNetworkDeviceConfiguration()
    net.attachment = VZNATNetworkDeviceAttachment()
    config.networkDevices = [net]

    // Input
    config.keyboards = [VZUSBKeyboardConfiguration()]
    config.pointingDevices = [VZUSBScreenCoordinatePointingDeviceConfiguration()]

    try config.validate()
    return config
}

// MARK: - App Delegate (GUI window via VZVirtualMachineView)

class AppDelegate: NSObject, NSApplicationDelegate {
    var window: NSWindow!
    var virtualMachine: VZVirtualMachine!
    var vmView: VZVirtualMachineView!

    func applicationDidFinishLaunching(_ notification: Notification) {
        do {
            let config = try createVMConfiguration()
            virtualMachine = VZVirtualMachine(configuration: config)

            vmView = VZVirtualMachineView()
            vmView.virtualMachine = virtualMachine
            vmView.automaticallyReconfiguresDisplay = true

            window = NSWindow(
                contentRect: NSRect(x: 0, y: 0, width: 1280, height: 800),
                styleMask: [.titled, .closable, .resizable, .miniaturizable],
                backing: .buffered,
                defer: false
            )
            window.title = "Ubuntu VM"
            window.contentView = vmView
            window.center()
            window.makeKeyAndOrderFront(nil)

            virtualMachine.start { result in
                switch result {
                case .success:
                    print("VM started. Close window to stop.")
                case .failure(let error):
                    fatalError("VM start failed: \(error)")
                }
            }
        } catch {
            fatalError("VM configuration failed: \(error)")
        }
    }

    func applicationShouldTerminateAfterLastWindowClosed(_ sender: NSApplication) -> Bool {
        true
    }
}

// MARK: - Main

let app = NSApplication.shared
app.setActivationPolicy(.regular)
let delegate = AppDelegate()
app.delegate = delegate
app.activate(ignoringOtherApps: true)
app.run()
