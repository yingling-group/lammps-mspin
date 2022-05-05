# Updating MSPIN

Instructions on how to keep the MSPIN source files updated to the lastest
stable version of LAMMPS.

## Get the Latest LAMMPS stable version files
To update to the latest LAMMPS version, clone the latest source code.

    git clone --depth 1 https://github.com/lammps/lammps.git --branch stable

    cd lammps

Set the path to the `lammps-mspin` repository as `GIT_WORK_TREE`.

    export GIT_WORK_TREE=/path/to/lammps-mspin

Check what has changed since the last update.
    
    git status

Overwrite the existing stable branch files.
    
    git reset --hard

## Update the MSPIN files
Open a new terminal session, make sure the `GIT_WORK_TREE` variable is clear.

Checkout the MSPIN files on top of the latest LAMMPS files
by repeating the steps as below.

    cd /path/to/lammps-mspin
    git status
    git reset --hard

Try building LAMMPS with the `mspin.cmake` preset. If successful, add and commit to git.
Otherwise, find the breaking changes, fix the conflicts and retry.
