#!/usr/bin/env python3
"""
Split VCF files into individual contact files and move them to the queb app
"""

import os
import re
import shutil
from pathlib import Path

def split_vcf_to_individual_contacts(vcf_file_path, output_dir):
    """Split a VCF file into individual contact files"""
    contacts = []
    current_contact = []
    
    try:
        with open(vcf_file_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except UnicodeDecodeError:
        with open(vcf_file_path, 'r', encoding='latin-1') as f:
            content = f.read()
    
    # Split by BEGIN:VCARD
    vcard_blocks = re.split(r'BEGIN:VCARD', content)
    
    for block in vcard_blocks[1:]:  # Skip first empty block
        if block.strip():
            contact_content = "BEGIN:VCARD" + block.strip()
            contacts.append(contact_content)
    
    # Create individual contact files
    base_name = os.path.splitext(os.path.basename(vcf_file_path))[0]
    contact_files = []
    
    for i, contact in enumerate(contacts):
        # Extract phone numbers for filename
        phone_match = re.search(r'TEL[^:]*:([^\n\r]+)', contact)
        if phone_match:
            phone = re.sub(r'[^\d+]', '', phone_match.group(1))
            if phone:
                filename = f"{base_name}_{phone}.vcf"
            else:
                filename = f"{base_name}_{i+1:03d}.vcf"
        else:
            filename = f"{base_name}_{i+1:03d}.vcf"
        
        # Ensure filename is safe
        filename = re.sub(r'[^\w\-_\.]', '_', filename)
        
        contact_file_path = os.path.join(output_dir, filename)
        
        with open(contact_file_path, 'w', encoding='utf-8') as f:
            f.write(contact)
        
        contact_files.append(contact_file_path)
    
    return contact_files

def main():
    # Create contacts directory in intention app
    contacts_dir = "/Users/m/Code/intention/contacts"
    os.makedirs(contacts_dir, exist_ok=True)
    
    # VCF files to process
    vcf_files = [
        "/Users/m/Code/Contacts.vcf",
        "/Users/m/Code/contacts_named.vcf", 
        "/Users/m/Code/contacts_no_name.vcf"
    ]
    
    all_contact_files = []
    
    for vcf_file in vcf_files:
        if os.path.exists(vcf_file):
            print(f"Processing {vcf_file}...")
            try:
                contact_files = split_vcf_to_individual_contacts(vcf_file, contacts_dir)
                all_contact_files.extend(contact_files)
                print(f"  Created {len(contact_files)} individual contact files")
            except Exception as e:
                print(f"  Error processing {vcf_file}: {e}")
        else:
            print(f"File not found: {vcf_file}")
    
    print(f"\nTotal contact files created: {len(all_contact_files)}")
    print(f"Contacts directory: {contacts_dir}")
    
    # List some example files
    if all_contact_files:
        print("\nExample contact files:")
        for i, file in enumerate(all_contact_files[:10]):
            print(f"  {os.path.basename(file)}")
        if len(all_contact_files) > 10:
            print(f"  ... and {len(all_contact_files) - 10} more")

if __name__ == "__main__":
    main()