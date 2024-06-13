<p align="left"><img src="https://rcedgar.github.io/usearch12_documentation/usearch12_banner.jpg" height="100"/></p>

Usearch implements several popular biological sequence search and clustering algorithms, including USEARCH, UCLUST, UPARSE, UCHIME, UNOISE and SINTAX.

Version 12 is the first open-source version of usearch. Compared to earlier versions, functionality which is well covered by other open-source projects has been removed. In particular, there is no support for OTU table manipulation or diversity analysis which is well supported by other tools such as [QIIME](https://qiime2.org/) and [DADA2](https://benjjneb.github.io/dada2/). The goal here is to simplify the package as much as reasonably possible to encourage collaborators to join the open-source project.

### Documentation

Docs. web site: [https://rcedgar.github.io/usearch12_documentation/](https://rcedgar.github.io/usearch12_documentation/)


Docs. source: [https://github.com/rcedgar/usearch12_documentation](https://github.com/rcedgar/usearch12_documentation)

## Installation

Download a binary (executable) file. There are no dependencies, so that is typically all you need to do. Make sure the execute bit is set if you are using Linux or OSX, and make sure the binary is in your PATH (or, you can type the full path name). For more details see [https://rcedgar.github.io/usearch12_documentation/install.html](https://rcedgar.github.io/usearch12_documentation/install.html).

## Building from source

### Windows

To build from the MSVC app, load the solution file `usearch12.sln` and select `Build` then `Rebuild Solution` from the main menu bar. 

To build from the command line, run `src/build_win.bash` from a command prompt. This requires that the MSVC build tools are in your PATH. The `build_win.bash` script (1) checks that there are no uncommitted changes to the repo, (2) generates a new `gitver.txt` file with the latest commit hash, and (3) runs `MSBuild.txt` to compile and link `usearch12.exe`.

### Linux

The primary development environment is Microsoft Visual C++ (MSVC). This means that the Linux `Makefile` is generated automatically from the MSVC project file `usearch12.vcxproj`. 

A pre-generated `Makefile` is included in the `src/` directory. This means that you can run `make` in the usual way. Generally, this `Makefile` should not be manually edited because changes will be lost the next time it is generated. 

Alternatively you can run the `src/build_linux.py` script, which (1) checks that there are no uncommitted changes to the repo, (2) generates a new `gitver.txt` file with the latest commit hash, (3) generates a new `Makefile` from `usearch12.vcxproj`, and (4) runs `MSBuild.exe` (MSVC's equivalent of `make`) to compile and link `usearch12.exe`.

### OSX

Currently building on OSX is not implemented. The Linux `Makefile` might work as-is, at most a few simple tweaks will be needed.
