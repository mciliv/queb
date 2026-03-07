#!/usr/bin/env python3

import os
import sys
import json
import uuid
import argparse
import subprocess
from pathlib import Path
from datetime import datetime

# Local storage for sessions
BASE = Path.home() / ".agent"
SESSIONS = BASE / "sessions"
BASE.mkdir(exist_ok=True)
SESSIONS.mkdir(exist_ok=True)

API_KEY = os.environ.get("OPENAI_API_KEY")
MODEL = "gpt-5"

def now():
    return datetime.utcnow().isoformat()

def new_session():
    sid = str(uuid.uuid4())
    path = SESSIONS / f"{sid}.json"
    data = {"id": sid, "created": now(), "messages": []}
    path.write_text(json.dumps(data, indent=2))
    return sid

def load_session(sid):
    path = SESSIONS / f"{sid}.json"
    if not path.exists():
        sys.exit(f"session not found: {sid}")
    return json.loads(path.read_text()), path

def save_session(path, data):
    path.write_text(json.dumps(data, indent=2))

def list_sessions():
    for f in sorted(SESSIONS.glob("*.json")):
        data = json.loads(f.read_text())
        print(f"{data['id']}  {data['created']}")

def read_files(paths):
    out = []
    for p in paths:
        path = Path(p)
        if path.is_file():
            try:
                out.append(f"\nFILE: {path}\n{path.read_text()}")
            except:
                pass
    return "\n".join(out)

def read_dirs(dirs):
    out = []
    for d in dirs:
        for path in Path(d).rglob("*"):
            if path.is_file():
                try:
                    out.append(f"\nFILE: {path}\n{path.read_text()}")
                except:
                    pass
    return "\n".join(out)

def call_responses_api(prompt, messages=None):
    """
    Uses the new OpenAI Responses API (/v1/responses)
    """
    import urllib.request

    # Convert chat messages to the format expected by 'input' in Responses API
    # if provided. Otherwise just use the string prompt.
    if messages:
        # Note: In the new API, 'input' can take the whole conversation
        payload = {
            "model": MODEL,
            "input": messages,
            "store": True
        }
    else:
        payload = {
            "model": MODEL,
            "input": prompt,
            "store": True
        }

    body = json.dumps(payload).encode()

    req = urllib.request.Request(
        "https://api.openai.com/v1/responses",
        data=body,
        headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {API_KEY}"
        }
    )

    try:
        with urllib.request.urlopen(req) as res:
            resp_data = json.loads(res.read())
            # Extract output_text or traverse the 'output' items
            # The Responses API returns a list of items in 'output'
            output_items = resp_data.get("output", [])
            text_content = ""
            for item in output_items:
                if item.get("type") == "message":
                    for content_part in item.get("content", []):
                        if content_part.get("type") == "output_text":
                            text_content += content_part.get("text", "")
            
            return text_content, resp_data
    except Exception as e:
        return f"Error calling API: {str(e)}", None

def extract_code(text):
    blocks = []
    lines = text.splitlines()
    inside = False
    buf = []
    for line in lines:
        if line.startswith("```"):
            if inside:
                blocks.append("\n".join(buf))
                buf = []
                inside = False
            else:
                inside = True
        elif inside:
            buf.append(line)
    return blocks

def run_code(code):
    try:
        p = subprocess.run(
            code,
            shell=True,
            capture_output=True,
            text=True,
            timeout=60
        )
        return p.stdout + p.stderr
    except Exception as e:
        return str(e)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--task")
    parser.add_argument("-f", "--files", nargs="*", default=[])
    parser.add_argument("-d", "--dirs", nargs="*", default=[])
    parser.add_argument("-i", "--id")
    parser.add_argument("-l", "--list", action="store_true")
    parser.add_argument("-x", "--exec", action="store_true")

    args = parser.parse_args()

    if args.list:
        list_sessions()
        return

    if args.id:
        data, path = load_session(args.id)
        sid = args.id
    else:
        sid = new_session()
        data, path = load_session(sid)

    context = ""

    if args.files:
        context += read_files(args.files)

    if args.dirs:
        context += read_dirs(args.dirs)

    if args.task:
        prompt = args.task + "\n" + context
        data["messages"].append({"role": "user", "content": prompt})

    if not data["messages"]:
        return

    print(f"[session {sid}]")

    # Use the new Responses API
    response_text, raw_response = call_responses_api(None, messages=data["messages"])

    print(response_text)

    # In the Responses API, we should ideally append the returned items to context
    # but for compatibility with the existing JSON session format:
    data["messages"].append({"role": "assistant", "content": response_text})
    save_session(path, data)

    if args.exec:
        codes = extract_code(response_text)
        for code in codes:
            print("\n[exec]")
            print(run_code(code))

if __name__ == "__main__":
    main()
