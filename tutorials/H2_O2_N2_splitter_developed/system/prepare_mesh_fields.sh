#!/bin/bash
set -euo pipefail

echo "Cleaning old mesh/sets..."
rm -rf constant/polyMesh
rm -rf 0/polyMesh
rm -rf 1e-* 2e-* 3e-* backup_0

echo "Generating base mesh..."
blockMesh


echo "backing up refineMeshDict and copying first refinement"
cp system/refineMeshDict system/refineMeshDict.bak

echo "Selecting refinement cells and refining"
cp system/refineMeshDict1 system/refineMeshDict
topoSet
refineMesh -overwrite

echo "Selecting refinement cells and refining"
cp system/refineMeshDict2 system/refineMeshDict
topoSet
refineMesh -overwrite

echo "Selecting refinement cells and refining"
cp system/refineMeshDict3 system/refineMeshDict
topoSet
refineMesh -overwrite

cp system/refineMeshDict.bak system/refineMeshDict
echo "Backing up initial fields..."
mkdir -p backup_0
cp -r 0/* backup_0/

echo "Applying initial fields..."
setFields

echo "Checking mesh..."
checkMesh

echo "Done."
