#! /bin/bash

# Place this in your MPL style lib directory.
# Then execute it from that directory and it will fetch the
# latest mysty sheets and figfun

names=( "" "-sub" "-12" "-quad" )

for i in "${names[@]}"
do
    echo mysty$i
    wget -O "mysty$i.mplstyle" https://raw.githubusercontent.com/saintsfan342000/MPL-a-al-SKK/master/mysty$i.mplstyle > /dev/null 2>&1
done

echo figfun.py
wget -O "../../../../figfun.py" https://raw.githubusercontent.com/saintsfan342000/MPL-a-al-SKK/master/figfun.py > /dev/null 2>&1


