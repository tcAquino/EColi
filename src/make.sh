#!/bin/bash
sed -i'.bak' "s/.*using namespace ecoli::model_.*/  using namespace ecoli::model_$1;/" "EColi.cpp"
make EColi
mv "../bin/EColi" "../bin/EColi_model_$1"
rm *.bak

