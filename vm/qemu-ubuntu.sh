#!/bin/zsh
# QEMU Ubuntu arm64 VM (macOS M-series HVF accel)
# Usage: ./qemu-ubuntu.sh {install|boot}
# Needs: brew install qemu

set -e

DISK="ubuntu-vm.qcow2"
ISO="$HOME/Downloads/ubuntu-25.10-desktop-arm64.iso"
FIRMWARE="/opt/homebrew/share/qemu/edk2-aarch64-code.fd"  # brew qemu

if [[ ! -f "$FIRMWARE" ]]; then
  echo "Install QEMU: brew install qemu"
  exit 1
fi

case "$1" in
  install)
    qemu-img create -f qcow2 $DISK 20G
    qemu-system-aarch64 \
      -M virt,accel=hvf \
      -cpu host \
      -smp 2 \
      -m 2G \
      -drive if=pflash,format=raw,readonly=on,file=$FIRMWARE \
      -drive file=$DISK,if=virtio \
      -cdrom "$ISO" \
      -device virtio-gpu-pci \
      -device qemu-xhci \
      -device usb-kbd \
      -device usb-tablet \
      -netdev user,id=net0 \
      -device virtio-net-pci,netdev=net0 \
      -display cocoa
    ;;
  boot)
    [[ -f $DISK ]] || { echo "Run 'install' first"; exit 1; }
    qemu-system-aarch64 \
      -M virt,accel=hvf \
      -cpu host \
      -smp 2 \
      -m 2G \
      -drive if=pflash,format=raw,readonly=on,file=$FIRMWARE \
      -drive file=$DISK,if=virtio \
      -device virtio-gpu-pci \
      -device qemu-xhci \
      -device usb-kbd \
      -device usb-tablet \
      -netdev user,id=net0 \
      -device virtio-net-pci,netdev=net0 \
      -display cocoa
    ;;
  *)
    echo "Usage: $0 {install|boot}"
    ;;
esac