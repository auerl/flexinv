#!/bin/bash

cd out
make
cd ..

cd inv
make
cd ..

cd mat/suwa
make 
cd ../..

cd ma/bowa
make 
cd ../..

cd  adp
make 
cd ..
