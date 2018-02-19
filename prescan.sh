######################################################################################################
#  ** Written by Jonatan Öström, jonatan.ostrom@gmail.com ** 
#  This script scanns fortran source directories and produces a file with dependencies to be included in a makefile. 
#  It is probably not very robust or portable. Works on Ubuntu 16. Probably only works when there is one module per file, but could work with any number. 
######################################################################################################

instructions="
<< INSTRUCTIONS >>
  
  This script creates the dependencies for a fortran makefile and writes then to 'prerequisites.makefile'
  The Makefile should contain the statement
  
srcdirs:=dir1 dir2
  
  somewhere, and
  
include prerequisites.makefile

  at the end. Running

./prescan.sh 
  
  creates the file 'prerequisites.makefile'. Alternatively the statement 'srcdirs:=dir1 dir2 ...' can be left out
  and dependency bewteen any sources in directories dir1, dir2, etc. can be created with:

./prescan.sh dir1 dir2
"

##########################################################################################################
# Input and input error handling
if [[ "$#" == "0" ]] # if no arguments use directories in makefile
    then 
    echo -e "\n-no arguments given, getting directories from makefile ..."
    source_places=$(grep -i '^srcdirs' Makefile |cut -d'=' -f 2)
    
    if [[ "$source_places" == "" ]]
    then 
        echo -e "\n << ERROR >> Directory list from makefile empty"
        echo -e " - No makefile?"
        echo -e " - No 'srcdirs:=dir1 dir2 dir3' statement?"
        exit
    fi
    echo -e "\n-directories from makefile: '$source_places', checking validity ... "
        
    for word in $source_places
    do
        if [[ -d $word ]]
        then 
            echo -e "\e[1m$word\e[0m is a directory"
        else 
            echo -e "\n\e[31mERROR\e[0m \e[1m$word\e[0m is \e[1mNOT\e[0m a directory"
            echo -e "\n$instrucitons"
            exit
        fi
    done
elif ! [[ -d "$1" ]] # if argument not a directory give instructions
then 
    echo -e "\nERROR Non-directory argument given."
    echo "$instructions"
    exit
else
    source_places=("$@") #to avoid command line arguments, write directories in the parethesis
fi

#exit
################################################################################################################################

#source_places=("$@") #to avoid command line arguments, write directories in the parethesis
ext=".f90"
filename="prerequisites.makefile"
bd=build

echo -e "\n-listing $ext files ..."

# create which[] array of all files for each place
i=0 
for place in ${source_places[*]}
do 
    which[$i]="$place/*$ext"
    i=$[$i+1]
done
echo -e "\n << \e[32mFILES\e[0m >>" ${which[*]}

echo -e "\n-scanning files for modules..."

# save the module names
i=0
for var in $(grep -i '^module' ${which[*]} | cut -d' ' -f 2 | cut -d! -f 1) #get word after "module" statement (file.f90:module hej!lol)
do 
    mod_arr[$i]=$var
    i=$[$i+1]
done
echo -e "\n << \e[32mMODULES\e[0m >>" ${mod_arr[*]}

echo -e "\n-creating source- and object-lists for found modules ..."

# save file name corresponding to modules 
i=0
for var in $(grep -i '^module' ${which[*]} | cut -d: -f 1 ) #get filename from grep output (file.f90:module hej!lol)
do 
    fil_arr[$i]=$var
    obj_arr[$i]="${var%.f90}.o"
    i=$[$i+1]
done

echo -e "\n << \e[32mFILES ORDERED\e[0m >>" ${fil_arr[*]}
echo -e "\n << \e[32mOBJECTS ORDERED\e[0m >>" ${obj_arr[*]}

length=$i


# start writing output file
echo -e "\n-creating '$filename' ..."
echo -e "\n### DEPENDENCIES______________________________\n\n" > $filename

echo -e "\n-checking dependencies and writing to '$filename' ..."
echo -e "\n << \e[32mDEPENDENCIES\e[0m >> "
i=0
while [[ $i -lt $length ]]
    do
    echo "$bd/$(basename ${obj_arr[$i]}):\\" >> $filename #write dependent object "astrid.o"
    for dep in $(grep -i '^ *use' ${fil_arr[$i]} | cut -d',' -f 1 | cut -d'!' -f 1 | grep -oi [use/].* | cut -d' ' -f 2) #list 'use'd modules
        do 
        j=0
        while [[ $j -lt $length ]]  
            do mod="${mod_arr[$j]}" #loop thgouht the module-list for each 'use'
            if [[ "$dep" == "$mod" ]] #if required module is found then:
                then echo "$bd/$(basename ${obj_arr[$j]})\\" >> $filename #write astrid.o's prerequisutes
                echo -e "\e[32m prerequisute\e[0m $dep\e[32m of\e[0m ${fil_arr[$i]}\e[32m found in\e[0m ${fil_arr[$j]}"
                break
            fi 
            j=$[$j+1] 
            
            if [[ $j == $length ]] # if required module is not found among sources, give an error
            then echo -e " << \e[31mERROR\e[0m >> used module '$dep' not in the sources, spelling misstake in dependent file '${fil_arr[$i]}' ?"
            echo "exiting"
            exit
            fi
        done
    done
    echo -e "\n\n"  >> $filename
    i=$[$i+1]
done 

echo "done"

######################################################################################################




