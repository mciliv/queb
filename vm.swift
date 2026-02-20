import Foundation
import Virtualization

func expandedURL(from path: String) -> URL {
    let expanded = (path as NSString).expandingTildeInPath
    return URL(fileURLWithPath: expanded)
}

// Set up global properties.
let virtualMachineConfiguration = VZVirtualMachineConfiguration()
virtualMachineConfiguration.cpuCount = 2
// 2 GB of memory. This is an arbitrary amount.
virtualMachineConfiguration.memorySize = 2 * 1024 * 1024 * 1024 

// Add a serial console and audio devices.
virtualMachineConfiguration.serialPorts = [createConsoleConfiguration()]


let outputStream = VZVirtioSoundDeviceOutputStreamConfiguration()
outputStream.sink = VZHostAudioOutputStreamSink()
let outputSoundDevice = VZVirtioSoundDeviceConfiguration()
outputSoundDevice.streams = [outputStream]


let inputStream = VZVirtioSoundDeviceInputStreamConfiguration()
inputStream.source = VZHostAudioInputStreamSource()
let inputSoundDevice = VZVirtioSoundDeviceConfiguration()
inputSoundDevice.streams = [inputStream]


virtualMachineConfiguration.audioDevices = [outputSoundDevice, inputSoundDevice]

// Paths for the Ubuntu ISO and a writable disk image.
// You can set an absolute path like "/Users/m/Downloads/ubuntu-25.10-desktop-arm64.iso"
// or a tilde path like "~/Downloads/ubuntu-25.10-desktop-arm64.iso".
let preferredISOPath = "~/Downloads/ubuntu-25.10-desktop-arm64.iso" // Change if needed.
var ubuntuISOURL = expandedURL(from: preferredISOPath)

// If the preferred path doesn't exist, auto-detect the first .iso in ~/Downloads.
if !FileManager.default.fileExists(atPath: ubuntuISOURL.path) {
    let downloadsURL = expandedURL(from: "~/Downloads")
    if let iso = try? FileManager.default.contentsOfDirectory(at: downloadsURL, includingPropertiesForKeys: nil)
        .first(where: { $0.pathExtension.lowercased() == "iso" }) {
        ubuntuISOURL = iso
    }
}

// Ensure the installer ISO exists.
guard FileManager.default.fileExists(atPath: ubuntuISOURL.path) else {
    fatalError("Ubuntu ISO not found. Set preferredISOPath correctly (e.g., /Users/yourname/Downloads/ubuntu-*.iso) or place an .iso in ~/Downloads.")
}

// Create (if needed) a 20 GB sparse disk image for installation target.
let supportDir = URL(fileURLWithPath: NSHomeDirectory()).appendingPathComponent("Library/Application Support/UbuntuVM", isDirectory: true)
try? FileManager.default.createDirectory(at: supportDir, withIntermediateDirectories: true)
let diskImageURL = supportDir.appendingPathComponent("disk.img")
if !FileManager.default.fileExists(atPath: diskImageURL.path) {
    // Create a sparse file of 20 GB.
    let fd = open(diskImageURL.path, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR)
    if fd == -1 {
        fatalError("Failed to create disk image at \(diskImageURL.path)")
    }
    // 20 * 1024 * 1024 * 1024 bytes
    let size: off_t = 20 * 1024 * 1024 * 1024
    if ftruncate(fd, size) != 0 {
        close(fd)
        fatalError("Failed to size disk image to 20GB")
    }
    close(fd)
}

// Attach the writable disk image as a Virtio block device.
let writableAttachment = try VZDiskImageStorageDeviceAttachment(url: diskImageURL, readOnly: false)
let writableDisk = VZVirtioBlockDeviceConfiguration(attachment: writableAttachment)

// Attach the Ubuntu ISO as a read-only CD-ROM device.
let isoAttachment = try VZDiskImageStorageDeviceAttachment(url: ubuntuISOURL, readOnly: true)
let isoDisk = VZVirtioBlockDeviceConfiguration(attachment: isoAttachment)

// Assign storage devices (writable disk first, then installer ISO).
virtualMachineConfiguration.storageDevices = [writableDisk, isoDisk]

// Add a Virtio graphics device so we can see the installer.
let graphicsDevice = VZVirtioGraphicsDeviceConfiguration()
let scanout = VZVirtioGraphicsScanoutConfiguration(widthInPixels: 1920, heightInPixels: 1080)
graphicsDevice.scanouts = [scanout]
virtualMachineConfiguration.graphicsDevices = [graphicsDevice]

// Validate this configuration against the current hardware;
// exit on error.
try! virtualMachineConfiguration.validate()

// Start the VM and check for successful startup when the block returns.
let virtualMachine = VZVirtualMachine(configuration: virtualMachineConfiguration)
virtualMachine.start(completionHandler: { (result) in
    switch result {
    case let .failure(error):
        fatalError("Virtual machine failed to start \(error)")
    default:
        NSLog("Virtual machine successfully started.")
    }
})

