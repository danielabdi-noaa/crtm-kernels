# Sample CRTM kernel

To comple with PGI compiler on hera, do this

    module load pgi/19.10 cuda/10.1

Then login to a compute node and compile and run as

    ./compile.sh
    ./runs.sh

Sample result

     > ./run.sh 
     Running version 1
     =======================
      Filling arrays with random values
      Finished filling arrays with random values.
      Creating RTV
      Finished creating RTV
      Calling kernel
      
      
      Finished executing kernel in =    3.635673    
      
        4.1067769403291986E-004   3.4815537764400261E-004   3.2094101497336910E-004 
        3.7737042782473498E-004   4.0001425534907497E-004   4.2665285055200813E-004 
        3.9698054551823953E-004   4.1921605504222445E-004   3.6573964432067582E-004 
        3.3485130584706490E-004   3.4071103146568200E-004   4.2450461947272412E-004 
        3.7793018510133558E-004   4.3977231668591708E-004   4.1620522915486143E-004 
        3.8670981495306835E-004
     Running version 2
     =======================
      Filling arrays with random values
      Finished filling arrays with random values.
      Creating RTV
      Finished creating RTV
      Calling kernel
      
      
      Finished executing kernel in =   2.4424000E-02
      
        4.1067769403291981E-004   3.4815537764400272E-004   3.2094101497336905E-004 
        3.7737042782473498E-004   4.0001425534907497E-004   4.2665285055200808E-004 
        3.9698054551823953E-004   4.1921605504222445E-004   3.6573964432067582E-004 
        3.3485130584706490E-004   3.4071103146568195E-004   4.2450461947272412E-004 
        3.7793018510133558E-004   4.3977231668591708E-004   4.1620522915486143E-004 
        3.8670981495306835E-004
     Running on cpu
     =======================
      Filling arrays with random values
      Finished filling arrays with random values.
      Creating RTV
      Finished creating RTV
      Calling kernel
      
      
      Finished executing kernel in =   0.4021940    
      
        4.1067769403291981E-004   3.4815537764400272E-004   3.2094101497336905E-004 
        3.7737042782473498E-004   4.0001425534907497E-004   4.2665285055200808E-004 
        3.9698054551823953E-004   4.1921605504222445E-004   3.6573964432067582E-004 
        3.3485130584706490E-004   3.4071103146568195E-004   4.2450461947272412E-004 
        3.7793018510133558E-004   4.3977231668591708E-004   4.1620522915486143E-004 
        3.8670981495306835E-004
