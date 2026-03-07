import Foundation
import Virtualization

func expandedURL(from path: String) -> URL {
    let expanded = (path as NSString).expandingTildeInPath
    return URL(fileURLWithPath: expanded)
}

let virtualMachineConfiguration = VZVirtualMachineConfiguration()
virtualMachineConfiguration.cpuCount = 2
virtualMachineConfiguration.memorySize = 2 * 1024 * 1024 * 1024 

let outputStream = VZVirtioSoundDeviceOutputStreamConfiguration()
outputStream.sink = VZHostAudioOutputStreamSink()
let outputSoundDevice = VZVirtioSoundDeviceConfiguration()
outputSoundDevice.streams = [outputStream]

let inputStream = VZVirtioSoundDeviceInputStreamConfiguration()
inputStream.source = VZHostAudioInputStreamSource()
let inputSoundDevice = VZVirtioSoundDeviceConfiguration()
inputSoundDevice.streams = [inputStream]

virtualMachineConfiguration.audioDevices = [outputSoundDevice, inputSoundDevice]

let supportDir = URL(fileURLWithPath: NSHomeDirectory()).appendingPathComponent("Library/Application Support/UbuntuVM", isDirectory: true)
try? FileManager.default.createDirectory(at: supportDir, withIntermediateDirectories: true)
let diskImageURL = supportDir.appendingPathComponent("disk.img")

guard FileManager.default.fileExists(atPath: diskImageURL.path) else {
    fatalError("Disk not found. Run vm.swift (installer) first to create/install.")
}

let writableAttachment = try VZDiskImageStorageDeviceAttachment(url: diskImageURL, readOnly: false)
let writableDisk = VZVirtioBlockDeviceConfiguration(attachment: writableAttachment)

virtualMachineConfiguration.storageDevices = [writableDisk]

let graphicsDevice = VZVirtioGraphicsDeviceConfiguration()
let scanout = VZVirtioGraphicsScanoutConfiguration(widthInPixels: 1920, heightInPixels: 1080)
graphicsDevice.scanouts = [scanout]
virtualMachineConfiguration.graphicsDevices = [graphicsDevice]

let networkDevice = VZVirtioNetworkDeviceConfiguration()
networkDevice.attachment = VZNATNetworkDeviceAttachment()
virtualMachineConfiguration.networkDevices = [networkDevice]

try! virtualMachineConfiguration.validate()

let virtualMachine = VZVirtualMachine(configuration: virtualMachineConfiguration)
virtualMachine.start { result in
  switch result {
  case .success:
    print("Ubuntu VM booted from disk. Ctrl+C to stop.")
  case .failure(let error):
    fatalError("VM start failed: \(error)")
  }
}

RunLoop.main.run()