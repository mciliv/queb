#!/usr/bin/env python3
import sys
import subprocess
import shlex
import re
import os

def run_gemini(prompt):
    subprocess.run(['gemini', '-p', prompt], text=True)

def run_search(term):
    subprocess.run(['ddgr', term])

def classify_input(user_input):
    prompt = f'''Classify ONLY as one line: \nSHELL:<exact cmd>\nLLM:<query>\nSEARCH:<term>\n\nInput: {user_input}'''
    try:
        result = subprocess.run(['gemini', '-p', prompt], capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            response = result.stdout.strip().split('\n')[0]
            return response
    except:
        pass
    return f'SHELL: {user_input}'

def execute_action(action, user_input):
    try:
        prefix, content = action.split(':', 1)
        content = content.strip()
    except ValueError:
        prefix, content = 'SHELL', user_input

    if prefix == 'SHELL':
        subprocess.run(shlex.split(content))
    elif prefix == 'LLM':
        run_gemini(content)
    elif prefix == 'SEARCH':
        run_search(content)
    else:
        subprocess.run(shlex.split(user_input))

def main():
    print("SmartCLI ready. Type 'exit' to quit.")
    while True:
        try:
            user_input = input('> ').strip()
            if user_input.lower() in ['exit', 'quit']:
                sys.exit(0)
            if not user_input:
                continue
            action = classify_input(user_input)
            print(f'Classified: {action}')
            execute_action(action, user_input)
        except EOFError:
            break
        except KeyboardInterrupt:
            print('\nGoodbye!')
            break
        except Exception as e:
            print(f'Error: {e}')

if __name__ == '__main__':
    main()
