breed [atoms atom]

atoms-own [
  ax     ; x-component of acceleration vector
  ay     ; y-component of acceeleration vector
  vx     ; x-component of velocity vector
  vy     ; y-component of velocity vector
]

globals [
  r-min
  eps
  -eps*4
  cutoff-dist
  force-offset
  PE-offset
  dt
  kb
]

;;;;;;;;;;;;;;;;;;;;;;
;; Setup Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;

to setup
  clear-all
  set-default-shape turtles "circle"
  set dt .01
  set kb (1 / 10)  ; just picking a random constant for Kb that makes things work reasonably
  set r-min  2 ^ (1 / 6)
  set eps 1
  set -eps*4 -4 * eps
  ifelse num-atoms = 3 [
    set cutoff-dist world-width
  ] [
    set cutoff-dist 2.5 * r-min
  ]
  set force-offset -1 * calc-force-without-offset cutoff-dist
  set PE-offset -1 * calc-pair-PE-without-offset cutoff-dist
  setup-atoms
  init-velocity
  if num-atoms = 3 [
    ask atoms [set vx 0 set vy 0]
    create-link-rulers
  ]

  reset-timer
  reset-ticks
end

to create-link-rulers
  ask atoms [
    ask other atoms [
      create-link-with myself
    ]
  ]
  ask links [set-label-distance]
end

to set-label-distance
  let d [distance [end1] of myself ] of end2
  set label (word precision (d / r-min) 1 " r0")
end

to setup-atoms
  create-atoms num-atoms [
    set shape "circle"
    set color blue
  ]

  if initial-config = "Solid" [
    let l sqrt(num-atoms) ;the # of atoms in a row
    let x-dist r-min
    let y-dist sqrt (x-dist ^ 2 - (x-dist / 2) ^ 2)
    let ypos (- l * x-dist / 2) ;the y position of the first atom
    let xpos (- l * x-dist / 2) ;the x position of the first atom
    let r-num 0  ;the row number
    ask turtles [  ;set the atoms; positions
      if xpos > (l * x-dist / 2)  [  ;condition to start a new row
        set r-num r-num + 1
        set xpos (- l * x-dist / 2) + (r-num mod 2) * x-dist / 2
        set ypos ypos + y-dist
      ]
      setxy xpos ypos  ;if we are still in the same row
      set xpos xpos + x-dist
    ]
  ]

  if initial-config = "Random" [
    ask atoms [
      setxy random-xcor random-ycor
    ]
    remove-overlap ;make sure atoms aren't overlapping
  ]

end

to init-velocity
  let v-avg sqrt( 2 * Kb * temp)
  ask atoms [
    let a random-float 1  ; a random amount of the total velocity to go the x direction
    set vx sqrt (a * v-avg ^ 2) * positive-or-negative
    set vy sqrt( v-avg ^ 2 - vx ^ 2)  * positive-or-negative
  ]
end


to-report positive-or-negative
  report ifelse-value random 2 = 0 [-1] [1]
end


to remove-overlap
  ask atoms [
    while [overlapping] [
      setxy random-xcor random-ycor
    ]
  ]
end

to-report overlapping
  report any? other turtles in-radius r-min
end


;;;;;;;;;;;;;;;;;;;;;;;;
;; Runtime Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;;;

to go
  (ifelse
    go-mode = "simulate" [simulate]
    go-mode = "drag atoms" [drag-atoms-with-mouse]
  )

end

to simulate
  ask links [hide-link]
  ask atoms [
    update-force-and-velocity
  ]
  ask atoms [
    move
  ]
  if constant-temp? [scale-velocities]
  tick-advance dt
  update-plots
end


to-report current-temp
  report (1 / (2 * Kb)) * mean [vx ^ 2 + vy ^ 2] of atoms
end


to scale-velocities
  let ctemp current-temp
  (ifelse
    ctemp = 0 and temp != 0 [
      ask atoms [init-velocity]
    ]
    ctemp != 0 [
      let scale-factor sqrt( temp / ctemp )  ; if "external" temperature is higher atoms will speed up and vice versa
      ask atoms [
        set vx vx * scale-factor
        set vy vy * scale-factor
      ]
    ]
  )
end

to update-force-and-velocity  ; atom procedure
  let new-ax 0
  let new-ay 0
  ask other atoms in-radius cutoff-dist [
    let r distance myself
    let force calc-force r
    face myself
    rt 180
    set new-ax new-ax + (force * dx)  ; assuming mass = 1. If not, neeed to divide force by mass to get acceleration
    set new-ay new-ay + (force * dy)
  ]
  set vx velocity-verlet-velocity vx ax new-ax
  set vy velocity-verlet-velocity vy ay new-ay
  set ax new-ax
  set ay new-ay

end

to-report calc-force [r]
  report calc-force-without-offset r + force-offset
end

to-report calc-force-without-offset [r]
  let r^3 r * r * r
  let r^6 r^3 * r^3
  ;report (-eps*4 / (r^6 * r)) * ((1 / r^6) - 1)
  report (-eps*4 * 6 / (r^6 * r)) * ((2 / r^6) - 1)
end

to move  ; atom procedure
  ;; Uses velocity-verlet algorithm
  set xcor velocity-verlet-pos xcor vx ax
  set ycor velocity-verlet-pos ycor vy ay
end

to-report velocity-verlet-pos [pos v a]  ; position, velocity and acceleration
  report pos + v * dt + (1 / 2) * a * (dt ^ 2)
end

to-report velocity-verlet-velocity [v a new-a]  ; velocity, previous acceleration, new acceleration
  report v + (1 / 2) * (new-a + a) * dt
end


to-report calc-PE
  let U 0  ;; U stands for PE

  ask other turtles in-radius cutoff-dist [

    set U U + calc-pair-PE (distance myself)
  ]
  report U
end

to-report calc-pair-PE [r]
  report calc-pair-PE-without-offset r + PE-offset
end

to-report calc-pair-PE-without-offset [r]
  let rsquare r ^ 2
  let attract-term 1 / rsquare ^ 3
  let repel-term attract-term * attract-term
  ;NOTE could do this a little faster by attract-term * (attract-term -1)
  report 4 * eps * (repel-term - attract-term)
end


to drag-atoms-with-mouse
  ask links [show-link]
  if mouse-down? [
    ask min-one-of turtles [distancexy mouse-xcor mouse-ycor] [
      let oldx xcor
      let oldy ycor
      setxy mouse-xcor mouse-ycor
      ifelse calc-PE > 1 [ ; if energy would be too high, don't let the atom go there.
        setxy oldx oldy
      ] [
        set vx 0
        set vy 0
      ]
    ]
  ]
  ask links [set-label-distance]
  display
end
@#$#@#$#@
GRAPHICS-WINDOW
180
10
514
345
-1
-1
29.69
1
18
1
1
1
0
1
1
1
-5
5
-5
5
1
1
1
ticks
30.0

BUTTON
105
55
170
90
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
105
105
172
138
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
0
10
172
43
num-atoms
num-atoms
0
30
20.0
1
1
NIL
HORIZONTAL

CHOOSER
0
50
95
95
initial-config
initial-config
"Solid" "Random"
1

SLIDER
0
155
172
188
temp
temp
0
8
0.3
.1
1
NIL
HORIZONTAL

SWITCH
0
195
135
228
constant-temp?
constant-temp?
0
1
-1000

MONITOR
0
235
55
280
temp
current-temp
3
1
11

CHOOSER
0
100
95
145
go-mode
go-mode
"drag atoms" "simulate"
1

TEXTBOX
60
235
170
275
This will match the temp slider if contant-temp? is on
11
0.0
1

@#$#@#$#@
## WHAT IS IT?

This is a molecular dynamics (MD) model using the Lennard-Jones potential function. This function models the fact that atoms attract each other when they are a small distance apart and repel each other when they are very close together. By modeling many atoms behaving according to the Lennard-Jones potential, we can see how the bulk behavior of matter at different temperatures emerges from the interactions between discrete atoms. The details of the Lennard-Jones function are discussed in the next section.


## HOW IT WORKS

MD simulations operate according to Newton's laws. The basic steps of the model are as follows. Each tick, each atom:

- Calculates the force that it feels from all other atoms using the Lennard-Jones potential
- Calculates its acceleration based on the net force and its mass using a = F / m
- Updates its velocity based on its acceleration
- Updates its position based on its velocity.

### The Lennard-Jones potential
The Lennard-Jones potential tells you the potential energy of an atom, given its distance from another atom. The derivative of the Lennard-Jones potential tells yout he force an atom feels from another atom based on their distance.

The potential is: V=4ϵ[(σ/r)^12−(σ/r)^6]. Where V is the intermolecular potential between two atoms or molecules, ϵ is depth of the potential well, σ is the distance at which the potential is zero (visualized as the diameter of the atoms), and r is the center-to-center distance of separation between both particles. This is an approximation; the potential function of a real atom depends on its electronic structure and will differ somewhat from the Lennard-Jones potential.
Atoms that are very close will be strongly pushed away from one another, and atoms that are far apart will be gently attracted. Make sure to check out the THINGS TO TRY section to explore the Lennard-Jones potential more in depth.

## HOW TO USE IT

### Simulation initialization
**initial-config**: select if you want the atoms to start out randomly positioned or in a hexagonally-close-packed structure.

**num-atoms**: select the number of atoms to start the simulation

**temp**: select the initial temperature (this will determine the average initial velocity of the atoms.

### Run the simulation

**constant-temp?**: turn this on if you want temperature to be held constant. This can be turned on and off during the simulation run. When it is on, you can move the **temp** slider to change the temperature.

## THINGS TO NOTICE

At low temperatures, the atoms will solidfy and at high temperatures they will break apart (evaporate). Also notice the the packing structure of atoms when they are in a solid and the structures that form when you cool the atoms down after evaporating compared to when they start in HCP.


## HOW TO CITE

If you mention this model or the NetLogo software in a publication, we ask that you include the citations below.

For the model itself:

* Kelter, J. and Wilensky, U. (2005).  NetLogo Lennard-Jones Molecular Dynamics model.  http://ccl.northwestern.edu/netlogo/models/Electrostatics.  Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

Please cite the NetLogo software as:

* Wilensky, U. (1999). NetLogo. http://ccl.northwestern.edu/netlogo/. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.


## COPYRIGHT AND LICENSE

Copyright 2021 Jacob Kelter and Uri Wilensky.

![CC BY-NC-SA 3.0](http://ccl.northwestern.edu/images/creativecommons/byncsa.png)

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 License.  To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 559 Nathan Abbott Way, Stanford, California 94305, USA.

Commercial licenses are also available. To inquire about commercial licenses, please contact Uri Wilensky at uri@northwestern.edu.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.3.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
1
@#$#@#$#@