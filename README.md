# optimize-temperature-set

## Setup Guide for Julia on Jupyter Notebook
This guide explains how to set up an environment to use the Julia language with Jupyter Notebooks (.ipynb files). It's ideal for running Julia code interactively in an editor like VSCode.

## Prerequisites
- **Julia**: Please install the version appropriate for your OS from the official website.
- **(Recommended) VSCode + Julia Extension**: This provides the smoothest development experience.

## Setup Steps:
### 1. Clone:
```zsh
git clone https://github.com/tatsuya0926/optimize-temperature-set.git
```

### 2. Initialize the Julia Project
Start the Julia interactive environment (REPL) by typing julia in your terminal or command prompt.

Type `julia` to launch Julia.<br>
Once Julia starts, press the `]` key to switch to the package mode, indicated by the `pkg>` prompt.
```Julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.x.x
 _/ |\__'_|_|_|\__'_|  |
|__/                   |

julia> ]

```
A project named after the current working directory will launch.
```Julia
(v1.x) pkg> activate .
  Activating project at `~/optimize-temperature-set`

(optimize-temperature-set) pkg>
```
### 3. Install Dependencies
Install all packages listed in the Project.toml file (such as `IJulia` and `Revise`).
```Julia
(optimize-temperature-set) pkg> instantiate
```
When `IJulia` is installed, a necessary Python environment (Miniconda) for Jupyter integration will be automatically downloaded and configured. **You generally do not need to install Python yourself.**<br>
Press the `Backspace` (or `Delete`) key to exit package mode and return to the `julia>` prompt.

### 4. Start Using Jupyter Notebook
The setup is now complete. You can start using it in two ways.
#### A) Using with VSCode (Recommended)
1. Open your project folder (optimize-temperature-set) in VSCode.
2. Choose a file with an `.ipynb` extension.
3. In the top right, click "Select Kernel" and choose your installed Julia version.
4. You can now write and execute Julia code in the cells.

#### B) Using with a Web Browser
Execute the following commands in the Julia REPL to launch the familiar Jupyter Notebook interface in your browser.
```Julia
julia> using IJulia
julia> notebook()
```

### (Appendix) Using an Existing Python Environment
If you already have a Python environment on your PC, like Anaconda, and wish to use it, you can specify which Python `PyCall` should use with the following steps.

1. Find the Python Path:<br>
In your terminal, run `which python` (macOS/Linux) or `where python` (Windows) and copy the full path to the Python executable you want to use.
2. Set the Path in Julia:<br>
In the Julia REPL, assign the copied path to `ENV["PYTHON"]`.
```Julia
# Example for macOS/Linux
julia> ENV["PYTHON"] = "/Users/yourname/anaconda3/bin/python"

# Example for Windows
julia> ENV["PYTHON"] = "C:\\Users\\yourname\\Anaconda3\\python.exe"
```
3. Rebuild PyCall:<br>
Enter package mode and rebuild `PyCall` to apply the new setting.
```
(optimize-temperature-set) pkg> build PyCall
```