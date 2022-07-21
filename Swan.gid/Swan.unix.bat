#!/bin/sh -f

rm -f "$2/$1.post.res"

# OutputFile: $2/$1.log
# ErrorFile: $2/$1.err

# delete the line before and uncomment the following line 
# to execute the program
mv "$2/$1.dat" "$2/$1.m"
KERNEL=`uname -s`
if [ $KERNEL != "Darwin" ]
then
#  matlab "$2/$1.m"
fi

