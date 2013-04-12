#sum_all.awk
{   for (i=1; i<=NF; i++) { sum[i]+= $i }   }

END { for (i=1; i<=NF; i++ ) 
	{ print i, sum[i] } 
          }
          
{   for (i=1; i<=NF; i++) { s += $i }   }

END { for (i=1; i<=1; i++ ) 
	{ print "Col", all, " =", s } 
          }