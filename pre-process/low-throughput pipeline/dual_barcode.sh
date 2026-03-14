#!/bin/bash
input=$1
# input contents:
# column1: OUTER_barcode
# column2: start INNER_barcode [1]
# column3: end INNER_barcode [96]
# column4: cell name prefix (optional)

pipdir=./
POOL=$pipdir/nano_96_barcode.txt
seq1=ATCT
seq2=CTACACGACGCTCTTCCGATCT

if [ -s barcode.fa ];then
	echo "Existed barcode.fa was removed!"
	rm barcode.fa
fi

while read i
do
    a=($i)
    if [ ${a[0]} -ge 0 ];then
        OUTER_BC=${a[0]}
        INNER_start=${a[1]}
        INNER_end=${a[2]}
        OUTER_Name=${a[3]}
        [ $INNER_start ] || INNER_start=1
        [ $INNER_end ] || INNER_end=96
        [ $OUTER_Name ] || OUTER_Name=O${OUTER_BC}
        OUTER_Seq=`cat $POOL | awk -v i=${OUTER_BC} '$1==i {print $2}'`
        for INNER_BC in `seq ${INNER_start} ${INNER_end}`
        do
            INNER_Seq=`cat $POOL | awk -v i=$INNER_BC '$1==i {print $2}'`
            if [ $INNER_BC -le 9 ];then
                Name=${OUTER_Name}I0${INNER_BC}
            else
                Name=${OUTER_Name}I${INNER_BC}
            fi
            BARCODE=${seq1}${OUTER_Seq}${seq2}${INNER_Seq}
            echo "Name: $Name"
            echo "OUTER_Seq: $OUTER_Seq"
            echo "INNER_Seq: $INNER_Seq"
            echo "=============================================================="
            echo ">${Name}" >> barcode.fa
            echo "$BARCODE" >> barcode.fa
        done
    else
        echo "Barcode is no integer!"
        exit 1
    fi
done < $input
