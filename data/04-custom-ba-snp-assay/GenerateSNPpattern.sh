## This script will take a SNP pattern file of the following format
#inputfile
##TTTTTCTCTC
##AAAAGAAAAA
##TTTCTTTTTT
##TTTCTTTTTT
##TTTCTTTTTT
##TTTTCCCTTT
##Gives an output file with the position number of the variant nucleotide along with the variant nucleotide (SNP).
##outputfile
##TTTTTCTCTC C 6 8 10 
##AAAAGAAAAA G 5 
##TTTCTTTTTT C 4 
##TTTCTTTTTT C 4 
##TTTCTTTTTT C 4 
##TTTTCCCTTT C 5 6 7 

awk '
    {
        n = split( $0, a, "" );
        for( i = 1; i <= n; i++ )
        {
            count[a[i]]++;
            pos[a[i]] = sprintf( "%s%d ", pos[a[i]], i );
        }

        min = "";
        for( x in count )
        {
            if( match( x, "[ACGT]" ) && (min == "" || count[x] < count[min] ) )
                min = x;
        }

        print $0, min, pos[min];

        delete count;
        delete pos;
    }
'   inputfile > outputfile
