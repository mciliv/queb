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
  ubuntu-vm

SHORTCUTS (in app menu bar)
  Cmd+R     Restart VM (wipes disk/EFI, boots fresh)
  Cmd+Q     Quit

DETAILS
  CPU/RAM   4 cores, 4 GB
  Disk      20 GB sparse image at ~/Library/Application Support/UbuntuVM/disk.img
  EFI       ~/Library/Application Support/UbuntuVM/efi-variable-store
  ISO       Auto-detects ~/Downloads/ubuntu-25.10-desktop-arm64.iso
            or first .iso in ~/Downloads
  Network   NAT (appears as wired ethernet in Ubuntu, not WiFi)

BUILD
  swiftc vm/vm.swift -o ubuntu-vm -framework Virtualization -framework Cocoa
  codesign --force --sign - --entitlements vm/entitlements.plist ubuntu-vm
""")
}

if command == "-h" || command == "--help" || command == "help" {
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

// MARK: - App Delegate

class AppDelegate: NSObject, NSApplicationDelegate {
    var window: NSWindow!
    var virtualMachine: VZVirtualMachine!
    var vmView: VZVirtualMachineView!

    func applicationDidFinishLaunching(_ notification: Notification) {
        setupMenuBar()
        startVM()
    }

    func setupMenuBar() {
        let mainMenu = NSMenu()

        // App menu
        let appMenuItem = NSMenuItem()
        let appMenu = NSMenu()
        appMenu.addItem(withTitle: "Quit Ubuntu VM", action: #selector(NSApplication.terminate(_:)), keyEquivalent: "q")
        appMenuItem.submenu = appMenu
        mainMenu.addItem(appMenuItem)

        // VM menu
        let vmMenuItem = NSMenuItem()
        let vmMenu = NSMenu(title: "VM")
        vmMenu.addItem(withTitle: "Restart VM", action: #selector(restartVM), keyEquivalent: "r")
        vmMenuItem.submenu = vmMenu
        mainMenu.addItem(vmMenuItem)

        NSApp.mainMenu = mainMenu
    }

    func startVM() {
        do {
            let config = try createVMConfiguration()
            virtualMachine = VZVirtualMachine(configuration: config)

            vmView = VZVirtualMachineView()
            vmView.virtualMachine = virtualMachine
            vmView.automaticallyReconfiguresDisplay = true

            if window == nil {
                window = NSWindow(
                    contentRect: NSRect(x: 0, y: 0, width: 1280, height: 800),
                    styleMask: [.titled, .closable, .resizable, .miniaturizable],
                    backing: .buffered,
                    defer: false
                )
                window.title = "Ubuntu VM"
                window.center()
            }
            window.contentView = vmView
            window.makeKeyAndOrderFront(nil)

            virtualMachine.start { result in
                if case .failure(let error) = result {
                    fatalError("VM start failed: \(error)")
                }
            }
        } catch {
            fatalError("VM configuration failed: \(error)")
        }
    }

    @objc func restartVM() {
        // Stop current VM, wipe state, start fresh
        if virtualMachine.canStop {
            virtualMachine.stop { [weak self] error in
                if let error { NSLog("Stop error: \(error)") }
                DispatchQueue.main.async { self?.wipeAndRestart() }
            }
        } else {
            wipeAndRestart()
        }
    }

    func wipeAndRestart() {
        vmView.virtualMachine = nil
        virtualMachine = nil
        try? FileManager.default.removeItem(at: supportDir)
        print("Cleared VM state. Restarting...")
        startVM()
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
