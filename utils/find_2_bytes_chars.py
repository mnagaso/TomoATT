import os
import re

def find_kanji_chars_in_file(file_path):
    kanji_pattern = re.compile(r'[\u4E00-\u9FFF]')
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
        content = file.read()
        matches = kanji_pattern.findall(content)
        if matches:
            print(f"Found Kanji characters in {file_path}: {''.join(matches)}")

def scan_directory(directory):
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            # skip non-text files
            target_files = ['.cpp', '.h', '.hpp', '.c', '.cc', '.hh', '.cxx', '.hxx', '.py', '.txt', '.md', '.rst', '.ipynb', '.sh']
            if not file_path.endswith(tuple(target_files)):
                continue

            find_kanji_chars_in_file(file_path)

if __name__ == "__main__":
    #directory_to_scan = '.'  # Change this to the directory you want to scan
    list_target_dir = ['./src', './test', './include', './examples']

    for directory_to_scan in list_target_dir:
        scan_directory(directory_to_scan)