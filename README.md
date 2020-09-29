# Privacy Preserving Structure-from-Motion

This is the COLMAP version for the Privacy Preserving SfM project.
The code is heavily based on standard [COLMAP](www.github.com/colmap/colmap),
but we removed a lot of code (likely we missed something) that is not directly
necessary for the privacy preserving pipeline to reduce the complexity and make
the code easier to understand.

### About

In this project we explore privacy preserving methods for SfM by replacing the
standard keypoint positions with lines to hide dynamic image content from
a processing server. The work was presented at the European Conference on
Computer Vision 2020.

If you are using this code in a scientific project please cite
```
@inproceedings{Geppert2020ECCV,
  author    = {Marcel Geppert and
               Viktor Larsson and
               Pablo Speciale and
               Johannes L. Sch{\"{o}}nberger and
               Marc Pollefeys},
  title     = {Privacy Preserving Structure-from-Motion},
  booktitle = {European Conference on Computer Vision (ECCV)},
  year      = {2020},
}

@inproceedings{schoenberger2016sfm,
  author={Sch\"{o}nberger, Johannes Lutz and Frahm, Jan-Michael},
  title={Structure-from-Motion Revisited},
  booktitle={Conference on Computer Vision and Pattern Recognition (CVPR)},
  year={2016},
}
```

Please not that this is not the exact same code that we used to create the results in our paper.

### Running the pipeline

The simplest way to create a privacy preserving model is to add the necessary
additional information (namely gravity direction and camera calibration) to the
input image folder and then run the pipeline as in standard COLMAP in the GUI
(feature extracton, matching, and reconstruction).

Providing the additional information can be done by adding two extra text files
per image as in
```
<image_name>.jpg
<image_name>.jpg.camera_model.txt
<image_name>.jpg.gravity.txt
```

### File formats

#### Gravity direction
The gravity direction is expected as a space separated text file with a single
line, the gravity vector pointing downwards in camera coordinates, e.g.
```
0.022717 0.999539 -0.020155
```
Note that we are using the computer vision convention for the coordinate system
with the x-axis pointing to the right,  y-axis downwards and z-axis to the
front (from the viewer into the image plane).


#### Camera calibration

For the camera calibration we also expect a text file with a
single line, containing first the camera model name and then the parameters
as comma separated values, e.g.
```
OPENCV 2575.939766, 2608.293891, 1599.263118, 1257.125449, 0.141865, -0.465301, 0.000000, 0.000000
```
Please check the `camera_models.h` file for the supported camera models, their
identifiers and the expected parameters.

### Installation

As our code is based on COLMAP, the build process is mostly identical.

> We deliberately removed the install targets to avoid conflicts with a
> potentially installed standard COLMAP.

We tested the following steps for Ubuntu 18.04:

1. Make sure you have all dependencies installed
    ```shell script
    sudo apt install \
        cuda \
        git \
        cmake \
        build-essential \
        libboost-program-options-dev \
        libboost-filesystem-dev \
        libboost-graph-dev \
        libboost-regex-dev \
        libboost-system-dev \
        libboost-test-dev \
        libeigen3-dev \
        libsuitesparse-dev \
        libfreeimage-dev \
        libgoogle-glog-dev \
        libgflags-dev \
        libglew-dev \
        qtbase5-dev \
        libqt5opengl5-dev 
    ```
   
2. Clone and build ceres
    ```shell script
    sudo apt install libatlas-base-dev libsuitesparse-dev
    git clone https://ceres-solver.googlesource.com/ceres-solver
    cd ceres-solver
    git checkout $(git describe --tags) # Checkout the latest release
    mkdir build
    cd build
    cmake .. -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF
    make -j
    sudo make install
    ```

3. Clone and build the pipeline
    ```shell script
    git clone --recurse-submodule https://github.com/colmap/privacy_preserving_sfm
    cd privacy_preserving_sfm
    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo
    make -j8
    ```
4. Run the gui
    ```shell script
    src/exe/ppsfm gui
    ```

### Datasets

We provide a dataset as an example input on our [Project Page](http://www.cvg.ethz.ch/research/privacy-preserving-sfm).


### Problems

We tried to minimize the complexity of the code and test the pipeline
extensively, but there is still a high chance that we missed certain cases.
If you encounter crashes or other problems while running the pipeline,
please let us know by creating an Issue here on GitHub.

