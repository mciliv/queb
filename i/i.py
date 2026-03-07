#!/usr/bin/env python3

import os
import sys
import argparse
import json
from pathlib import Path

import pyperclip
import openai
openai.api_key = os.getenv("OPENAI_API_KEY")

def main():
    content = ' '.join(content) if (content := sys.argv[1:]) else pyperclip.paste()
    chat = Path("cure.json")
    message = {"role": "user", "content": content}

    completion = openai.ChatCompletion.create(
      model="gpt-4",
      messages=new_message(chat, message))

    response = completion.choices[0].message
    print(response.content)
    new_message(chat, response)

def new_message(chat, message) -> list[dict]:
    if chat.exists():
        with chat.open("r") as f:
            file_content = json.load(f)
        if isinstance(file_content, list):
            file_content.append(message)
    else:
        file_content = [message]
    with chat.open("w") as f:
        json.dump(file_content, f)
    return file_content

if __name__ == '__main__':
    main()

