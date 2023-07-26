# A Shape Modulus for Fractal Geometry Generation

## Intro

![Four intricate fractal shapes in the shape of an armadillo, bunny, dragon, and human hand](https://github.com/alexaschor/JuliaShapeModulus/assets/77604978/9946f6d4-413a-427e-a83d-78d05163340b)

This is the code for the paper titled _A Shape Modulus for Fractal Geometry Generation_, published at SGP2023 by Alexa Schor (me) and [Theodore Kim](https://tkim.graphics).

Download the paper: [PDF (192MiB)](http://tkim.graphics/VERSOR/SchorKim_2023.pdf), [PDF compressed (4MiB)](http://tkim.graphics/VERSOR/SchorKim_2023_smaller.pdf)


_Abstract:_ We present an efficient new method for computing Mandelbrot-like fractals (Julia sets) that approximate a user-defined shape. Our algorithm is orders of magnitude faster than previous methods, as it entirely sidesteps the need for a time-consuming numerical optimization. It is also more robust, succeeding on shapes where previous approaches failed. The key to our approach is a versor-modulus analysis of fractals that allows us to formulate a novel shape modulus function that directly controls the broad shape of a Julia set, while keeping fine-grained fractal details intact. Our formulation contains flexible artistic controls that allow users to seamlessly add fractal detail to desired spatial regions, while transitioning back to the original shape in others. No previous approach allows Mandelbrot-like details to be "painted" onto meshes.

## This Repo

First of all, thanks for coming to this page! If you have any questions or ideas about the paper or code or any problems with building, running or modifying the code, please email me at <img height="20px" alt="an image of an email address. the emaill address is the same as my github username, followed by at gmail.com" src="https://github.com/alexaschor/JuliaShapeModulus/assets/77604978/1225837c-eda3-4873-b92b-3cebe3f42d2e">
 and I'll do my best to help. 

This is the code used to generate the meshes used in most of the paper figures. For an interactive, (mostly) real time 2D version that allows you to play with the parameters, check out https://github.com/alexaschor/JSM2D

If you want to run Kim 2015's method, the [original code release](http://www.tkim.graphics/JULIA/source.html) is now difficult to build because it relies on the old C++ version of BEM++ (now called bempp), which is now maintained in python. I have a [modified version](https://github.com/alexaschor/QUIJIBO) hosted on this GitHub account that should be identical, just ported over to the new version and with a compiler error fixed.

## Compilation
It should build with `make`, though you might have to edit `projects/include.mk` with your compiler and compiler flags. This code has no external dependencies with the exception of one optional one: there's a parallel mesh compute script (symlinked in `bin/`) that runs N separate jobs to compute different parts of the mesh and then uses trimesh2's `mesh_cat` to join them together.

## Usage
Compiling the project produces three executables, compiled separately from directories in `projects/`. They each have their own Makefile, but each of them includes `projects/include.mk`, so that's where you can set the compiler to use and any flags you want to specify. 

The three executables spit out documentation when called with no arguments, which is reproduced here for convenience.

### 1. `genPoly`, for generating polynomial versor functions
```
❯ ./bin/genPoly
USAGE:
-----------------------------------------------
To inspect a *.p4d polynomial:
 ./bin/genPoly INS <*.p4d file>
-----------------------------------------------
To generate a random polynomial in a bounding box:
 ./bin/genPoly GEN <numRoots> <min power> <max power> <*.p4d output> BOX <x min> <y min> <z min> <x max> <y max> <z max>
To generate a random polynomial in an OBJ's bounding box:
 ./bin/genPoly GEN <numRoots> <min power> <max power> <*.p4d output> OBJBOX <*.obj model>
To generate a random polynomial near an SDF's surface:
 ./bin/genPoly GEN <numRoots> <min power> <max power> <*.p4d output> SDF <*.f3d sdf> <max distance>
-----------------------------------------------
To generate a polynomial whose roots are the union of two input polynomials:
 ./bin/genPoly ADD <first *.p4d input> <second *.p4d input> <*.p4d output>
-----------------------------------------------
To scale the positions of a polynomial's roots in space:
 ./bin/genPoly SCALE <*.p4d input> <factor> <*.p4d output>
```

### 2. `sdfGen`, for generating SDFs from a triangle mesh (essentially the same as [christopherbatty/SDFGen](https://github.com/christopherbatty/SDFGen))
```
❯ ./bin/sdfGen
USAGE:
To generate an SDF from a mesh with automatically generated bounds:
 ./bin/sdfGen <*.obj input> <resolution> <*.f3d output> <padding cells>
To get the bounds for a mesh sequence:
 ./bin/sdfGen BOUNDS <obj 1> <obj 2> ... <obj N>
To generate an SDF from a mesh with specified bounds:
 ./bin/sdfGen <*.obj input> <resolution> <*.f3d output> <min X> <min Y> <min Z> <max X> <max Y> <max Z>
```

### 3. `run`, for computing the shape-modulus Julia set given a polynomial for the versor function and an SDF
A note on this one: if you're just trying to reproduce the results from the paper, you probably don't need to use any offset for the origin. It's an interesting parameter to play with and causes some cool effects, but unless you want to tweak it just set the offset to `0 0 0`. The octree specifier string is for if you want to compute just a subsection of the overall Julia set, it's useful for zooming in on something or for running in parallel. If you want to run things in parallel, check out `bin/prun`. It's a symlink to a Python script that will automatically run this in parallel, and if you have Trimesh2's `mesh_cat` in your PATH it will also join the resultant meshes back together.
```
USAGE:
To create a shaped julia set from a distance field:
 ./bin/run <SDF *.f3d> <P(q) *.poly4d, or the word RANDOM> <output resolution> <a> <b> <offset x> <offset y> <offset z> <output *.obj> <optional: octree specifier string>

    This will compute the 3D Julia quaternion Julia set of the function:
        f(q) = r(q) * d(q)
    Where:
        r(q) = e^( a * SDF(q) + b )
        d(q) = Normalized( P(q) )

    P(q) is the quaternion polynomial function, which determines the character of the fractal detail.
    a is a parameter which controls the thickness of the shell in which the chaotic effect has significant influence
    on set membership, and b is a parameter which controls the position along the surface of the SDF of that shell.

    The offset X, Y, and Z parameters move the origin of the dynamical system around in space, which causes
    the Julia set to dissolve in interesting ways.
    The octree specifier string is an optional parameter useful for computing large Julia sets in parallel.
    You can select a box in an evenly-subdivided octree of arbitrary depth specified by a string of digits 0-7,
    laid out as follows:
                 +---+
               / |4|5|
              /  +-+-+
             /   |7|6|
            /    +---+
           +---+    /
           |0|1|   /
           +-+-+  /
           |3|2| /
        ▲  +---+
        |
        Y X--▶
        Z ●
        (into page)
    Each character of the string will go one level deeper, so the string '5555' specifies the 1/16-edge length box at the far back corner.
```
## Usage (as a library)

Aside from the main building executables, this repository can also be used as a library for representing vector and scalar functions on R3 and computing Julia sets in general. The main strength of the library is the ability to represent scalar and vector fields in the same way whether they are defined by some implicit function and evaluated on-the-fly or they are stored in-memory. Caching, interpolation, gradient calculation, file I/O, and meshing via Marching Cubes with root-finding are all supported. The scalar and vector field code are in `src/field.h`, and should be pretty much self-contained. The Julia set-specific ones are in `src/julia.h`, and are self-contained except for the QUATERNION and POLYNOMIAL-4D classes in `lib/Quaternion`. 
