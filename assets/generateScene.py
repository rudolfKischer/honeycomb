

#sample sphere
# {
#     "name": "rightWall",
#     "position": [502.0, 0.0, -4.0],
#     "radius": 500,
#     "color": [1.0, 0.0, 0.0],
#     "emission": [0.0, 0.0, 0.0],
#     "orm": [0.0, 1.0, 0.0]
# },

# orm: oppacity, roughness, metallic

import json;
import os;
from pathlib import Path

def roughnessMetalGrid(n):

  color = [1.0, 1.0, 1.0]

  scene = {}
  scene["spheres"] = []
  for i in range(n):
    for j in range(n):
      sphere = {}
      sphere["name"] = "sphere" + str(i) + str(j)
      sphere["radius"] = 2
      sphere["position"] = [i*sphere["radius"]*2 - (n * sphere["radius"]) / 2, j*sphere["radius"]*2 - (n * sphere["radius"]) / 2, -4.0]
      sphere["color"] = color
      sphere["emission"] = [0.0, 0.0, 0.0]
      sphere["orm"] = [0.0, (n-i)/n, j/n]
      scene["spheres"].append(sphere)
  
  return scene

def opacityRoughnessGrid(n):

  color = [1.0, 1.0, 1.0]

  scene = {}
  scene["spheres"] = []
  for i in range(1,n):
    for j in range(1,n):
      sphere = {}
      sphere["name"] = "sphere" + str(i) + str(j)
      sphere["radius"] = 2
      sphere["position"] = [i*sphere["radius"]*2 - (n * sphere["radius"]), j*sphere["radius"]*2 - (n * sphere["radius"]), -4.0]
      sphere["color"] = [(i)/n, (n-i)/n, j/n] #color
      sphere["emission"] = [0.0, 0.0, 0.0]
      sphere["orm"] = [0.0,1.0,0.0]#[i/n, j/n, 0.0]
      scene["spheres"].append(sphere)

  return scene

def cornellBox():
  # use really large spheres to create walls, floor and ceiling
  # make the ceiling emit light
  wallR = 1000
  scene = {}
  scene["spheres"] = []
  room_length = 10
  directions = [[1.0, 0.0, 0.0], 
                [-1.0, 0.0, 0.0], 
                [0.0, 1.0, 0.0], 
                [0.0, -1.0, 0.0], 
                [0.0, 0.0, 1.0], 
                [0.0, 0.0, -1.0]]
  for i in range(6):
    sphere = {}
    sphere["name"] = "wall" + str(i)
    sphere["radius"] = wallR
    sphere["position"] = [ wallR * v for v in directions[i]]
    room = [ room_length * v for v in directions[i]]
    sphere["position"] = [ room[j] + sphere["position"][j] for j in range(3)]
    sphere["color"] = [1.0, 1.0, 1.0]
    sphere["emission"] = [0.0, 0.0, 0.0]
    sphere["orm"] = [0.0, 1.0, 0.0]
    scene["spheres"].append(sphere)


  # make right wall red
  scene["spheres"][0]["color"] = [1.0, 0.0, 0.0]

  # make left wall green
  scene["spheres"][1]["color"] = [0.0, 1.0, 0.0]

  # make the back wall emissive 
  scene["spheres"][4]["emission"] = [0.2, 0.2, 0.2]

  # add a light in the ceiling
  sphere = {}
  sphere["name"] = "ceiling"
  sphere["radius"] = 20.0;
  sphere["position"] = [0.0, room_length - 1.0 + sphere["radius"], 0.0]
  sphere["color"] = [1.0, 1.0, 1.0]
  sphere["emission"] = [1.0, 0.9, 0.9]
  sphere["orm"] = [0.0, 1.0, 0.0]

  scene["spheres"].append(sphere)


  return scene





def main():
  # sceneGrid = roughnessMetalGrid(4)
  sceneGrid = opacityRoughnessGrid(4)
  sceneBox = cornellBox()

  # add the grid scene spheres to the cornell box
  scene = sceneBox
  for sphere in sceneGrid["spheres"]:
    scene["spheres"].append(sphere)


  localDir = Path(__file__).parent
  # file_path = f'{localDir}/grid.json'
  file_path = f'{localDir}/grid.json'

  with open(file_path, 'w') as outfile:
    json.dump(scene, outfile)

if __name__ == "__main__":
  main()


