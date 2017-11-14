# 3DCONF, relor and resect

Some few edits of the programs 3DCONF, relor and resect written by Bon A. Dewitt in 1998-2000.

## What ?

### 3DCONF

This program performs a three-dimensional conformal coordinate transformation by least squares.

### relor

Relative orientation program.

### resect

Single photo space resection program.

## Why ?

Minor change to ensure unix compatibility.
The purpose is to stop using Windows/DOS libraries `malloc.h` and `conio.h` that are only used for UI on Windows.

If you are working on windows you can still use the source codes in the *original* directory.


## Building

### Building with gcc

Command to build **3DCONF** :

    gcc -Wall 3DCONF.c -o 3DCONF

Command to build **relor** :

    gcc -Wall relor.c -o relor

Command to build **resect** :

    gcc -Wall resect.c -o resect

### Building with clang

Command to build **3DCONF** :

    clang 3DCONF.c -o 3DCONF

Command to build **relor** :

    clang relor.c -o relor

Command to build **resect** :

    clang resect.c -o resect

## How to use ?

The descriptions and instructions are defined in the comments at the end of each source code.

## Preprocessor constants

In these softwares, the maximum amount of points that can be processed are defined by preprocessor constants:

* In **3DCONF** :
    * `MAXCOM` at line 20
    * `MAXUNK` at line 21
* In **relor** :
    * `MAXPTS` at line 20
* In **resect** :
    * `MAXPTS` at line 19
