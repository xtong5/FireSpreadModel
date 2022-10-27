dirs="`ls | grep "frames_" | tr '\n' ' ' | cat`"
rootDir="`pwd`"

# ensure all dims the same and get frame counts
frameCounts=""
dims=""
for dir in $dirs; do
  cd $dir
  dims="${dims} `identify \`find . -name "*.png" -print | head -n 1\` | grep -o "[0-9]*x[0-9]* "`"
  frameCounts="${frameCounts} `ls | wc -l`"
  cd ..
done

maxFrames=`echo $frameCounts | tr ' ' '\n' | sort -n | tail -1`

dim="`echo $dims | tr ' ' '\n' | head -n 1 | tr -d '\n'`"
allSame=true
for thisDim in $dims; do
  if [[ $thisDim != $dim ]]; then
    allSame=false
  fi
done

if $allSame; then
width=`echo $dim | grep -o "[0-9]*" | head -n 1`
height=`echo $dim | grep -o "[0-9]*" | tail -n 1`

for dir in $dirs; do
  cd $dir
  ffmpeg -i frame_fireLine_%04d.png -vf scale=$width:-2 fire.mp4
  cd ..
done

else
echo "Not all dimensions agreed"
echo $dims

fi

python3 "4x4 vid.py" $width $height
cd frames
ffmpeg -i fire4x4_%04d.png -vf scale=$(( $width*2 )):-2 fire4x4.mp4
cd ..
