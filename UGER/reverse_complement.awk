 BEGIN{ 
     j= n = split("A C G T", t); 
     for (i=0; ++i <= n;) 
         map[t[i]] = t[j--]; 
     map["N"]="N";
 }
 NR%4==2{
     for (i=length; i; i--) 
	 printf "%s",map[substr($0,i,1)]; 
    print x 
 }
 NR%4==0{
     for (i=length; i; i--) 
	 printf "%s",substr($0,i,1); 
    print x 
 } 
 NR%4==1 ||  NR%4==3{
    print
 }
