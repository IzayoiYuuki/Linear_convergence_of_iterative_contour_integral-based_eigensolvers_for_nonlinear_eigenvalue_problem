# Convergence Theory of NLFEAST
This respository includes codes of the journal paper 'Linear convergence of iterative contour integral-based eigensolvers for nonlinear eigenvalue problems'.


All programs in this repository **must** be run with `./ConvNLFEAST-master` as the root directory. Before running any program, please make sure that:

1. You are currently in the `./ConvNLFEAST-master` directory.
2. All folders under this directory have been temporarily added to the MATLAB path.

To add the folders to the MATLAB path:

1. **If you are using the MATLAB graphical interface**:  
   Right-click `ConvNLFEAST-master` > `Add to Path` > `Selected Folders and Subfolders`.

2. **If you are using the command line**:  
   Open MATLAB in the `./ConvNLFEAST-master` directory and run

   ```matlab
   addpath(genpath('./'));
   ```

To test with photonics problem, you also need the MATLAB Symbolic Math Toolbox.
```matlab
https://ww2.mathworks.cn/products/symbolic.html
```

## Content
1. `./figure/` contains outputs, both pdf and txt.
1. `./lib/` contains external library files such as NLEVP.
1. `./matrix/` contains some test matrices.
1. `./paper/` contains generators of the figrues in the paper.
1. `./source/` contains eigensolvers such as NLFEAST and Beyn.

## How to reproduce the figures in the paper
1. Figure 1:
    1. run `./paper/test_intro.m` to get data.
    1. run `./paper/plot_intro.m` for plotting.
1. Table 2 and Fugure 3:
    1. run `./paper/9examples/test_all.m` to get data **(Take long time)**
    1. run `./paper/9examples/plot_timetable.m` to summarize the time consumed
    1. run `./paper/9examples/plot_all.m` to plot the convergence history.
1. Figure 4: 
    1. run `./paper/test_time2accuracy.m` to get data.
    1. run `./paper/plot_time2accuracy.m` for plotting.
1. Table 3 and Figure 5:
    1. run `./paper/same_eigv/test_same_eigv.m` to get data.
    1. run `./paper/same_eigv/plot_same_eigv.m` to get Figure 5.
1. Figure 6:
    1. run `./paper/NLFEAST_relaxed/test_NLFEAST_relaxed.m` to get data.
    1. run `./paper/NLFEAST_relaxed/plot_NLFEAST_relaxed.m` for plotting. 
    
