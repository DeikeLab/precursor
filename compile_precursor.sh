#
mkdir -p precursor
#
# compile
#
CC99='mpicc -std=c99' qcc -grid=octree -events -autolink -O2 -Wall -pipe -D_FORTIFY_SOURCE=2 -D_MPI=1 -DTREE -o precursor/precursor precursor.c -L$BASILISK/gl -lfb_tiny -lm
