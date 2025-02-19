import os
import re
import shutil
import subprocess

def replace_in_file(file_path, old_str, new_str):
    """
    Function to replace a string in the specified file.
    """
    try:
        # Open the file in read mode and read its content
        with open(file_path, 'r', encoding='utf-8') as file:
            content = file.read()

        # Perform the replacement
        if old_str in content:
            content = content.replace(old_str, new_str)
            # Open the file in write mode and save the updated content
            with open(file_path, 'w', encoding='utf-8') as file:
                file.write(content)
            print(f"Processed: {file_path}")
        else:
            print(f"No match found in: {file_path}")
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")

def process_directory(directory, file_extension, old_str, new_str):
    """
    Function to recursively search for files in the specified directory and perform replacements.
    """
    if not os.path.isdir(directory):
        print(f"Error: Directory '{directory}' does not exist.")
        return

    # Recursively search for all files in the directory
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(file_extension):
                file_path = os.path.join(root, file)
                replace_in_file(file_path, old_str, new_str)

def copy_directory(src_dir, dst_dir):
    """
    Recursively copy a directory and its contents to a new location.
    """
    # Check if the destination directory already exists
    if os.path.exists(dst_dir):
        print(f"Destination directory '{dst_dir}' already exists. Skipping copy.")
        return

    if not os.path.isdir(src_dir):
        print(f"Error: Source directory '{src_dir}' does not exist.")
        return
    try:
        # Copy the directory and its contents recursively
        shutil.copytree(src_dir, dst_dir)
        print(f"Successfully copied '{src_dir}' to '{dst_dir}'.")
    except FileExistsError:
        print(f"Error: Destination directory '{dst_dir}' already exists.")
    except Exception as e:
        print(f"An error occurred: {e}")

def replace_relative_paths(file_path):
    """
    Replace '../' with '../../' in the specified file.
    If '../' appears consecutively, only add one additional level.
    """
    try:
        # Read the file content
        with open(file_path, 'r', encoding='utf-8') as file:
            content = file.read()

        # Use regex to find and replace '../'
        # Pattern explanation:
        # - `(../)+`: Matches one or more consecutive '../'
        # - Replace with the same number of '../' plus one additional '../'
        def replace_match(match):
            # Count the number of '../' in the match
            count = len(match.group(0)) // 3  # Each '../' is 3 characters long
            return '../' * (count + 1)  # Add one more '../'

        updated_content = re.sub(r'(../)+', replace_match, content)

        # Write the updated content back to the file
        with open(file_path, 'w', encoding='utf-8') as file:
            file.write(updated_content)

        print(f"Processed: {file_path}")
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")


def main():
    src_directory = "../../free_energy/mbar_analysis/"
    target_directory = "../mbar_analysis"
    # Copy the mbar_analysis source code
    copy_directory(src_directory, target_directory)

    # Settings for the target directory and replacement strings
    file_extension = ".fpp"
    old_string = "ma_"
    new_string = "mbar_"

    # Execute the processing
    process_directory(target_directory, file_extension, old_string, new_string)

    # create Makefile.depends
    current_dir = os.getcwd()
    os.chdir(target_directory)
    command = ["make", "depend"]
    result = subprocess.run(command, check=True, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(result.stdout)
    print(result.stderr)
    os.chdir(current_dir)
    print("All .fpp files have been processed.")

    # replace_relative_paths(target_directory + "/Makefile")

if __name__ == "__main__":
    main()
