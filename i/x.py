import os
import argparse

from xai_sdk import Client
from xai_sdk.chat import user, system

parser = argparse.ArgumentParser()
parser.add_argument('system')
parser.add_argument('user')
parser.add_argument('--no-stream', action='store_true')
args = parser.parse_args()

client = Client(
    api_key=os.getenv("XAI_API_KEY"),
    timeout=3600, # Override default timeout with longer timeout for reasoning models
)

chat = client.chat.create(model="grok-4")
chat.append(system(args.system))
chat.append(user(args.user))

if args.no_stream:
    response = chat.sample()
    print(response.content)
else:
    for response, chunk in chat.stream():
        print(chunk.content, end="", flush=True) # Each chunk's content
        print(response.content, end="", flush=True) # The response object auto-accumulates the chunks

