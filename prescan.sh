######################################################################################################
#  ** Written by Jonatan Öström, jonatan.ostrom@gmail.com ** 
#  This script scanns fortran source directories and produces a file with dependencies to be included in a makefile. 
#  It is probably not very robust or portable. Works on Ubuntu 16. Probably only works when there is one module per file. Maybe no module in a file is ok.  
######################################################################################################

instructions="  Usage:

./prescan.sh src/ src2/ 
  
  Makefile dependencies from fortran sources in given directories are written to  'prerequisites.makefile' 
  Include it in the makefile with the statement 'include prerequisites.makefile' "
# If run with anything else than a directory argument, instructions are given. 
if ! [[ -d "$1" ]]
then echo "$instructions"
exit
fi

######################################################################################################

source_places=("$@") #to avoid command line arguments, write directories in the parethesis
ext=".f90"
filename="prerequisites.makefile"
bd=build

echo "listing $ext files ..."

# create which[] array of all files for each place
i=0 
for place in ${source_places[*]}
do which[$i]="$place/*$ext"
i=$[$i+1]
done

echo "scanning files ..."

# save the module names
i=0
for var in $(grep -i '^module' ${which[*]} | cut -d' ' -f 2 | cut -d! -f 1) 
do 
mod_arr[$i]=$var
i=$[$i+1]
done


# save file name corresponding to modules 
i=0
for var in $(grep -i '^module' ${which[*]} | cut -d: -f 1 ) 
do 
fil_arr[$i]=$var
obj_arr[$i]="${var%.f90}.o"
i=$[$i+1]
done

length=$i


objects=""
for dep in "${mod_arr[*]}" 
    do 
    j=0
    while [[ $j -lt $length ]]  
        do mod="${mod_arr[$j]}"
        if [[ "$dep" == "$mod" ]]
        then
            objects="$objects ${obj_arr[$j]}" 
        fi
        j=$[$j+1]
    done
done

# start writing output file
echo "creating '$filename' ..."
echo -e "\n# DEPENDENCIES______________________________\n" > $filename

i=0
while [[ $i -lt $length ]]
    do
    echo "$bd/$(basename ${obj_arr[$i]}):\\" >> $filename
    for dep in $(grep -i '^ *use' ${fil_arr[$i]} | cut -d',' -f 1 | grep -oi [use/].* | cut -d' ' -f 2)
        do 
        j=0
        while [[ $j -lt $length ]]  
            do mod="${mod_arr[$j]}"
            if [[ "$dep" == "$mod" ]]
                then echo "$bd/$(basename ${obj_arr[$j]})\\" >> $filename
            fi
            j=$[$j+1]
        done
    done
    echo -e "\n\n"  >> $filename
    i=$[$i+1]
done 

echo "done"

######################################################################################################




