import numpy as np
from PIL import Image, ImageDraw, ImageFont
from os import listdir, scandir, makedirs
from sys import argv

# init
i = 1
width, height = (3130, 1493)
if len(argv) == 3:
    width = int(argv[1])
    height = int(argv[2])
print(width, height, argv)
frame = np.empty((height*2, width*2, 3), np.uint8)
unfinished = [True] * 4
font = ImageFont.truetype('/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf', 72)
frames = np.empty((4, height, width, 3), np.uint8)

# find what values were used
subdirs = [dir.name for dir in scandir('.') if dir.is_dir() and 'frames_' in dir.name]
subscripts = [dirname.split('_')[1] for dirname in subdirs]

numFrames = max([len(listdir(f'{dir}')) for dir in subdirs])

makedirs('frames', exist_ok=True)

while any(unfinished):
    # extract frames
    for j, dir in enumerate(subdirs):
        try: 
            with Image.open(f'{dir}/frame_fireLine_{i:04}.png') as img:
                frames[j] = np.array(img)
        except FileNotFoundError: unfinished[j] = False

    # add frames into single 4x4 frame
    frame[0:height, 0:width, :] = frames[0]
    frame[0:height, width:, :] = frames[1]
    frame[height:, 0:width, :] = frames[2]
    frame[height:, width:, :] = frames[3]

    # convert array into PIL image
    img = Image.fromarray(frame)
    # annotate 4x4 with property value
    draw = ImageDraw.Draw(img)
    draw.text((0, 0), f"strengthSlopeDiffuse = {subscripts[0]}", (0,0,0), font=font)
    draw.text((width, 0), f"strengthSlopeDiffuse = {subscripts[1]}", (0,0,0), font=font)
    draw.text((0, height), f"strengthSlopeDiffuse = {subscripts[2]}", (0,0,0), font=font)
    draw.text((width, height), f"strengthSlopeDiffuse = {subscripts[3]}", (0,0,0), font=font)

    # save image
    img.save(f'frames/fire4x4_{i:04d}.png')
    img.close()
    print(100 * i // numFrames, '%', sep='')
    i += 1
