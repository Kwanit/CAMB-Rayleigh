"""
Python script to run and update Jupyter notebooks (CAMBdemo.ipynb or ScalEqs.ipynb).
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path


def check_command_exists(command):
    """Check if a command exists in PATH."""
    try:
        subprocess.run([command, "--version"], capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def check_python_package(package):
    """Check if a Python package is installed."""
    try:
        subprocess.run([sys.executable, "-c", f"import {package}"], capture_output=True, check=True)
        return True
    except subprocess.CalledProcessError:
        return False


def install_package(package):
    """Install a Python package using pip."""
    print(f"Installing {package}...")
    try:
        subprocess.run([sys.executable, "-m", "pip", "install", package], check=True)
        return True
    except subprocess.CalledProcessError:
        print(f"Failed to install {package}")
        return False


def execute_notebook(notebook_path, timeout=600):
    """Execute a Jupyter notebook and return the path to the executed notebook."""
    notebook_path = Path(notebook_path)
    temp_output_path = notebook_path.parent / f"{notebook_path.stem}_executed{notebook_path.suffix}"

    print(f"Running {notebook_path}...")

    try:
        cmd = [
            "jupyter",
            "nbconvert",
            "--to",
            "notebook",
            "--execute",
            str(notebook_path),
            "--output-dir",
            str(notebook_path.parent),
            "--output",
            temp_output_path.name,
            f"--ExecutePreprocessor.timeout={timeout}",
        ]
        subprocess.run(cmd, check=True)
        return temp_output_path
    except subprocess.CalledProcessError as e:
        print(f"Error executing notebook: {e}")
        return None


def check_notebook_errors(notebook_path):
    """Check if the executed notebook has any errors."""
    try:
        with open(notebook_path, encoding="utf-8") as f:
            notebook = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error reading notebook: {e}")
        return True

    for cell in notebook.get("cells", []):
        if cell.get("cell_type") == "code" and cell.get("outputs"):
            for output in cell["outputs"]:
                if output.get("output_type") == "error":
                    print("Error found in notebook execution:")
                    if "traceback" in output:
                        for line in output["traceback"]:
                            print(line)
                    return True
    return False


def clean_notebook(notebook_path):
    """Clean the notebook by removing trailing semicolons."""
    scripts_dir = Path(__file__).parent
    semicolon_script = scripts_dir / "remove_trailing_semicolons.py"

    if semicolon_script.exists():
        print("Cleaning up the notebook...")
        try:
            subprocess.run([sys.executable, str(semicolon_script), str(notebook_path)], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Warning: Failed to clean notebook: {e}")
    else:
        print("Warning: remove_trailing_semicolons.py script not found. Skipping semicolon cleanup.")


def generate_html(notebook_path):
    """Generate HTML version of the notebook."""
    notebook_path = Path(notebook_path)
    print("Generating HTML version of the notebook...")

    try:
        subprocess.run(
            ["jupyter", "nbconvert", "--to", "html", str(notebook_path), "--output-dir", str(notebook_path.parent)],
            check=True,
        )
        print("HTML version generated successfully.")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Warning: Failed to generate HTML version: {e}")
        return False


def main():
    """Main function to process notebook."""
    parser = argparse.ArgumentParser(description="Run and update Jupyter notebooks (CAMBdemo.ipynb or ScalEqs.ipynb)")
    parser.add_argument(
        "notebook",
        nargs="?",
        default="CAMBdemo",
        help="Notebook to process: 'CAMBdemo', 'ScalEqs', or path to notebook file (default: CAMBdemo)",
    )
    parser.add_argument("--timeout", type=int, default=600, help="Execution timeout in seconds (default: 600)")
    parser.add_argument("--no-html", action="store_true", help="Skip HTML generation")

    args = parser.parse_args()

    # Determine notebook path
    if args.notebook.endswith(".ipynb"):
        notebook_path = Path(args.notebook)
    elif args.notebook.lower() == "camb" or args.notebook.lower() == "cambdemo":
        notebook_path = Path("docs/CAMBdemo.ipynb")
    elif args.notebook.lower() == "scal" or args.notebook.lower() == "scaleqs":
        notebook_path = Path("docs/ScalEqs.ipynb")
    else:
        # Try to find the notebook in docs directory
        notebook_path = Path(f"docs/{args.notebook}.ipynb")
        if not notebook_path.exists():
            notebook_path = Path(f"docs/{args.notebook}")

    if not notebook_path.exists():
        print(f"Error: Notebook {notebook_path} not found.")
        sys.exit(1)

    # Check dependencies
    if not check_command_exists("jupyter"):
        print("Error: Jupyter is not installed or not in PATH. Please install it with 'pip install jupyter'.")
        sys.exit(1)

    # Check nbconvert
    try:
        subprocess.run(["jupyter", "nbconvert", "--version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Jupyter nbconvert is not installed. Installing it now...")
        if not install_package("nbconvert"):
            print("Error: Failed to install nbconvert. Please install it manually with 'pip install nbconvert'.")
            sys.exit(1)

    # Create backup
    backup_path = notebook_path.parent / f"{notebook_path.stem}_backup{notebook_path.suffix}"
    print(f"Creating backup of {notebook_path} to {backup_path}...")
    shutil.copy2(notebook_path, backup_path)

    # Execute notebook
    temp_output_path = execute_notebook(notebook_path, args.timeout)
    if not temp_output_path or not temp_output_path.exists():
        print("Error: Failed to execute notebook. Output file not found.")
        sys.exit(1)

    # Check for errors
    if check_notebook_errors(temp_output_path):
        print("Error: Notebook execution had errors. Keeping original notebook.")
        temp_output_path.unlink()
        sys.exit(1)

    # Replace original with executed notebook
    print(f"Execution successful. Updating {notebook_path} with executed notebook...")
    shutil.move(str(temp_output_path), str(notebook_path))

    # Clean notebook
    clean_notebook(notebook_path)

    # Generate HTML
    if not args.no_html:
        generate_html(notebook_path)

    print(f"Done! {notebook_path.name} has been successfully executed and updated.")


if __name__ == "__main__":
    main()
