#!/usr/bin/env python3
"""
Apple Contacts vCard Processor
This script helps process vCard files from Apple Contacts and identify entries with "No Name"
"""

import os
import re
import sys
from pathlib import Path

def parse_vcard(file_path):
    """Parse a vCard file and return a list of contact dictionaries"""
    contacts = []
    current_contact = {}
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except UnicodeDecodeError:
        with open(file_path, 'r', encoding='latin-1') as f:
            content = f.read()
    
    # Split by BEGIN:VCARD
    vcard_blocks = re.split(r'BEGIN:VCARD', content)
    
    for block in vcard_blocks[1:]:  # Skip first empty block
        contact = {}
        lines = block.strip().split('\n')
        
        for line in lines:
            if line.startswith('FN:') or line.startswith('N:'):
                # FN is formatted name, N is structured name
                key = 'FN' if line.startswith('FN:') else 'N'
                value = line[3:].strip()
                contact[key] = value
            elif line.startswith('TEL:'):
                if 'TEL' not in contact:
                    contact['TEL'] = []
                contact['TEL'].append(line[4:].strip())
            elif line.startswith('EMAIL:'):
                if 'EMAIL' not in contact:
                    contact['EMAIL'] = []
                contact['EMAIL'].append(line[6:].strip())
            elif line.startswith('ORG:'):
                contact['ORG'] = line[4:].strip()
            elif line.startswith('NOTE:'):
                contact['NOTE'] = line[5:].strip()
        
        if contact:  # Only add non-empty contacts
            contacts.append(contact)
    
    return contacts

def identify_no_name_contacts(contacts):
    """Identify contacts that have no name or 'No Name'"""
    no_name_contacts = []
    named_contacts = []
    
    for contact in contacts:
        # Check if contact has no name or is explicitly "No Name"
        has_name = False
        name_value = ""
        
        if 'FN' in contact and contact['FN'].strip():
            name_value = contact['FN'].strip()
            if name_value.lower() not in ['no name', 'noname', '']:
                has_name = True
        elif 'N' in contact and contact['N'].strip():
            name_value = contact['N'].strip()
            if name_value.lower() not in ['no name', 'noname', '']:
                has_name = True
        
        if not has_name or name_value.lower() in ['no name', 'noname', '']:
            no_name_contacts.append(contact)
        else:
            named_contacts.append(contact)
    
    return no_name_contacts, named_contacts

def write_vcard(contacts, output_file):
    """Write contacts to a vCard file"""
    with open(output_file, 'w', encoding='utf-8') as f:
        for contact in contacts:
            f.write("BEGIN:VCARD\n")
            f.write("VERSION:3.0\n")
            
            if 'FN' in contact:
                f.write(f"FN:{contact['FN']}\n")
            if 'N' in contact:
                f.write(f"N:{contact['N']}\n")
            if 'ORG' in contact:
                f.write(f"ORG:{contact['ORG']}\n")
            if 'TEL' in contact:
                for tel in contact['TEL']:
                    f.write(f"TEL:{tel}\n")
            if 'EMAIL' in contact:
                for email in contact['EMAIL']:
                    f.write(f"EMAIL:{email}\n")
            if 'NOTE' in contact:
                f.write(f"NOTE:{contact['NOTE']}\n")
            
            f.write("END:VCARD\n")

def main():
    print("Apple Contacts vCard Processor")
    print("=" * 40)
    
    # Look for vCard files in common locations
    possible_locations = [
        "/Users/m/Library/Containers/com.apple.AddressBook/Data/tmp/TemporaryItems/SharingItems/",
        "/Users/m/Desktop/",
        "/Users/m/Downloads/",
        "/Users/m/",
    ]
    
    vcard_files = []
    for location in possible_locations:
        if os.path.exists(location):
            for file in os.listdir(location):
                if file.endswith('.vcf') or file.endswith('.vcard'):
                    vcard_files.append(os.path.join(location, file))
    
    if not vcard_files:
        print("No vCard files found in common locations.")
        print("\nTo export your contacts from Apple Contacts:")
        print("1. Open the Contacts app")
        print("2. Select all contacts (Cmd+A)")
        print("3. Go to File > Export > Export vCard...")
        print("4. Save the file to your Desktop or Downloads folder")
        print("5. Run this script again")
        return
    
    print(f"Found {len(vcard_files)} vCard file(s):")
    for i, file in enumerate(vcard_files):
        size = os.path.getsize(file)
        print(f"  {i+1}. {file} ({size} bytes)")
    
    # Process each vCard file
    for vcard_file in vcard_files:
        if os.path.getsize(vcard_file) == 0:
            print(f"\nSkipping {vcard_file} (empty file)")
            continue
            
        print(f"\nProcessing: {vcard_file}")
        try:
            contacts = parse_vcard(vcard_file)
            print(f"Found {len(contacts)} contacts")
            
            no_name_contacts, named_contacts = identify_no_name_contacts(contacts)
            
            print(f"  - Named contacts: {len(named_contacts)}")
            print(f"  - No Name contacts: {len(no_name_contacts)}")
            
            if no_name_contacts:
                # Create output files
                base_name = os.path.splitext(os.path.basename(vcard_file))[0]
                no_name_file = f"{base_name}_no_name.vcf"
                named_file = f"{base_name}_named.vcf"
                
                write_vcard(no_name_contacts, no_name_file)
                write_vcard(named_contacts, named_file)
                
                print(f"\n  Created separate files:")
                print(f"    - {no_name_file} ({len(no_name_contacts)} contacts)")
                print(f"    - {named_file} ({len(named_contacts)} contacts)")
                
                # Show details of no-name contacts
                print(f"\n  No Name contacts details:")
                for i, contact in enumerate(no_name_contacts[:10]):  # Show first 10
                    print(f"    {i+1}. FN: {contact.get('FN', 'N/A')}, N: {contact.get('N', 'N/A')}")
                    if 'TEL' in contact:
                        print(f"       Phone: {', '.join(contact['TEL'])}")
                    if 'EMAIL' in contact:
                        print(f"       Email: {', '.join(contact['EMAIL'])}")
                
                if len(no_name_contacts) > 10:
                    print(f"    ... and {len(no_name_contacts) - 10} more")
            
        except Exception as e:
            print(f"Error processing {vcard_file}: {e}")

if __name__ == "__main__":
    main()