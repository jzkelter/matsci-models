globals[MAX-MOVE-DIST CUTOFF-DIST Vtotal EPS TOTAL-MOVE-ATTEMPTS TOTAL-SUCCESSFUL-MOVES DIAMETER POT-OFFSET CURRENT-MOVE-ATTEMPTS CURRENT-SUCCESSFUL-MOVES]

to setup
  clear-all
  reset-ticks
  resize-world -16 16 -16 16
  set EPS 1
  ;Set the diameter of particles based on density
  set DIAMETER sqrt(density * world-width * world-height / N) ;This is number density
  crt N [
    set shape "circle" 
    set size DIAMETER 
    set color blue 
  ]
  setup-atoms
  set MAX-MOVE-DIST DIAMETER
  set CUTOFF-DIST 2.5 * DIAMETER 
  set POT-OFFSET (- (4 * ((DIAMETER / CUTOFF-DIST) ^ 12 - (DIAMETER / CUTOFF-DIST) ^ 6)) )
  set Vtotal calc-Vtotal
  tick ;so that the tune-acceptance rate doesn't get called on first tick
  
end

to go

  ;Each tick, attempt N moves. On average, every particle moves each tick
  repeat N [  
    ask one-of turtles[
      attempt-move
    ]
  ]
  
  ;tune the move distance to adjust the acceptance rate every N ticks
  if ticks mod N = 0[
    tune-acceptance-rate
  ]
  
  tick
end


to attempt-move
  set TOTAL-MOVE-ATTEMPTS TOTAL-MOVE-ATTEMPTS + 1  ;the is the total running average
  set CURRENT-MOVE-ATTEMPTS CURRENT-MOVE-ATTEMPTS + 1 ;this is just since the last max-move-distance adjustment
  let Vold calc-V; Calculate current energy
  let deltax (random-float 2 * MAX-MOVE-DIST) - MAX-MOVE-DIST  ; pick random x distance
  let deltay (random-float 2 * MAX-MOVE-DIST) - MAX-MOVE-DIST ; pick random y distance
  setxy (xcor + deltax) (ycor + deltay) ;move the random x and y distances
  let Vnew calc-V ;Calcualte the new energy
  
  let deltaV Vnew - Vold
  ifelse (Vnew < Vold) or (random-float 1 < exp( - deltaV / T) )
    [
      set TOTAL-SUCCESSFUL-MOVES TOTAL-SUCCESSFUL-MOVES + 1   ;the is the total running average
      set CURRENT-SUCCESSFUL-MOVES CURRENT-SUCCESSFUL-MOVES + 1   ;this is just since the last max-move-distance adjustment
      set Vtotal Vtotal + deltaV
    ]
    [
      setxy (xcor - deltax) (ycor - deltay) ;reset position
    ] 
end


to-report calc-Vtotal
  let Vtot 0
  ask turtles [
    set Vtot Vtot + calc-V
  ]
  report Vtot / 2 ;divide by two because each particle has been counted twice
end

to-report calc-V
  let V 0
  let Vi 0
  let repel-term 0
  let attract-term 0
  let rsquare 0
  ask other turtles in-radius CUTOFF-DIST [
    ;set rsquare (xcor - [xcor] of myself) ^ 2 + (ycor - [ycor] of myself) ^ 2
    set rsquare (distance myself) ^ 2
    let dsquare DIAMETER * DIAMETER   
    set attract-term (dsquare * dsquare * dsquare)  / (rsquare * rsquare * rsquare)
    set repel-term attract-term * attract-term
    ;NOTE could do this a little faster by attract-term * (attract-term -1)
    set Vi 4 * EPS * (repel-term - attract-term) + POT-OFFSET
    set V V + Vi
  ]
  report V 
end

to-report accept-rate
  report CURRENT-SUCCESSFUL-MOVES / CURRENT-MOVE-ATTEMPTS 
end
 
to tune-acceptance-rate
  ifelse accept-rate < 0.5[
    set MAX-MOVE-DIST MAX-MOVE-DIST * .95
  ][
  set MAX-MOVE-DIST MAX-MOVE-DIST * 1.05
  if MAX-MOVE-DIST > DIAMETER [
    set MAX-MOVE-DIST DIAMETER
  ]
  ]
  set CURRENT-SUCCESSFUL-MOVES 0
  set CURRENT-MOVE-ATTEMPTS 0
  
end


to-report energy-per-particle
  report Vtotal / N
end
;*********SETUP Procedures*************

to setup-atoms
  if initial-config = "hcp"[
    let L sqrt(N) 
    let row-dist (2 ^ (1 / 6)) * DIAMETER ;this is the distance with minimum energy
    let ypos (- L * row-dist / 2)
    let xpos (- L * row-dist / 2)
    let rnum 0
    ask turtles [
      if xpos > (L * row-dist / 2)  [
        set rnum rnum + 1
        set xpos (- L * row-dist / 2) + (rnum mod 2) * row-dist / 2
        set ypos ypos + row-dist
      ]
      setxy xpos ypos
      set xpos xpos + row-dist
    ]
  ]

  if initial-config = "random"[
    ask turtles [
      setxy random-xcor random-ycor 
    ]
  ]
  
  remove-overlap
end

to remove-overlap
  let RMIN  0.7 * DIAMETER
  ask turtles [
    while [overlapping RMIN] [
      setxy random-xcor random-ycor 
    ]
    
  ]
end


to-report overlapping [RMIN]
  ifelse any? other turtles in-radius RMIN
  [report True]
  [report False]
end

@#$#@#$#@
GRAPHICS-WINDOW
240
10
724
515
16
16
14.364
1
10
1
1
1
0
1
1
1
-16
16
-16
16
1
1
1
ticks
30.0

BUTTON
30
195
204
228
NIL
setup\n
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
120
331
204
364
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
30
72
202
105
N
N
1
1000
497
1
1
NIL
HORIZONTAL

SLIDER
31
294
205
327
T
T
.01
2
0.3
.01
1
NIL
HORIZONTAL

SLIDER
30
109
202
142
density
density
0
.6
0.26
.01
1
NIL
HORIZONTAL

CHOOSER
30
146
203
191
Initial-Config
Initial-Config
"hcp" "random"
0

PLOT
33
368
215
513
energy per particle
NIL
NIL
0.0
10.0
-2.0
2.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "if ticks > 3 [\nplot vtotal / N\n]"

TEXTBOX
4
10
24
46
1
30
15.0
1

TEXTBOX
31
10
238
73
Model starting point. You can choose the number of atoms, the density and the intial configuration (random or hexagonally-close-packed)
11
0.0
1

TEXTBOX
31
229
208
257
_____________________________
11
0.0
1

BUTTON
32
331
117
364
go-once
go\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
2
251
20
287
2
30
15.0
1

TEXTBOX
30
248
225
290
Adjust the temparture and run the model. The temperature can be adjusted while the model runs\n
11
0.0
1

@#$#@#$#@
## WHAT IS IT?

Matter usually exists in one of three phases: solid, liquid or gas. We often think of matter as a continuous bulk substance that somewhat mysteriously changes between these phases as temperature changes (and as pressure changes). However, matter is actually made of discrete elements called atoms. The behavior of matter at different tempartures can be understood by understanding how atoms interact with each other. 

This model shows how the bulk behavior of matter at different temperatures emerges from the interactions between discrete atoms. 

## HOW IT WORKS

The basic principle of the model is simple: each atom tries to minimize its energy. Atoms are modeled as soft spheres (this means they can overlap). The energy of an atom depends soley on its distance from other atoms. The amount of energy one atom gives another is determined by a model potential function called the Lennard-Jones potential. The potential function of a real atom depends on its electronic structure and will differ somewhat from the Lennard-Jones potential. The following is a brief description of how atoms behave given the Lennard-Jones potential. 

The Lennard-Jones potential is: V=4ϵ[(σ/r)^12−(σ/r)^6]. Where V is the intermolecular potential between two atoms or molecules, ϵ is depth of the potential well, σ is the distance at which the potential is zero (visualized as the diameter of the atoms), and r is distance of separation between both particles (measued from their centers).

Atoms that are too close (r<σ) give each other positive energy, meaning that they will tend to move away from each other. After a certain distance (r>σ), atoms give each other negative energy. The energy of two atoms is minimized at a critical distance. 

The steps for each tick in the model are:

1. Choose a random atom.
2. Ask that atom to calculate its present energy. The atom does this by summing up the Lennard-Jones potentials with all the atoms within CUTTOFF-DISTANCE of itself (in nature all atoms interact with all others, but beyond a certain distance the effect is so negligble we can ignore it).
3. The atom then moves a random distance. It moves in both the x and y directions a distance less than or equal to its diameter.
4. The atom calculates its energy in the new tentative position.
5a. Finally, the atom compares its new tentative energy with the old energy. If the new energy is lower, it remains in the new position. 
5b. If the new energy is higher, it still may stay in the new position. The probability of staying in the new position, even if the new energy is higher, depends on the Temperature. The higher the Temperature, the more likely the atom is to accept a new, higher energy position.
6. Repeat steps 1-5 N times. N is the number of atoms. This means that on average each atom attempts one move each tick, but on any given tick some atoms may attempt multiple moves, and some may not attempt any.

The *T* slider controls the temperature. The higher the temperature the more likely an atom is to maintain a new position that is actually a higher energy configuration.

Note that the number of atoms is very small compared to most real systems. Also, real systems are three-dimensional, while this model is 2D.

To learn more about the Lennard-Jones potential see: 


* http://chemwiki.ucdavis.edu/Physical_Chemistry/Quantum_Mechanics/Atomic_Theory/Intermolecular_Forces/Lennard-Jones_Potential.
* http://en.wikipedia.org/wiki/Lennard-Jones_potential

## HOW TO USE IT

### 1) Simulation starting point

**number of atoms**: You can change the number of atoms by adjusting the **N** slider.

**density**: You can adjust the density by adjusting the **density** slider. This wil not effect the number of atoms. It is will change their size to change the density. Density is calculated as a number density (number per unit area) as opposed to area density (the summed area of all atoms per total area). 


**initial configuration**: You can choose the intial configuraiton of the atoms. They can either start randomly distributed, or in a low energy hexogonally close packed (HCP) structure.

### 2) Run the model
**temperature**: Set an initial temperature before running the model. You can change the temperature during the run as well. 

You can run the model one step at a time with the **go once** button. This will go through steps 1-6 above once. Or, you can have the model run continuously with the **go** button. 

## THINGS TO NOTICE

At high temeratures the atoms move around the environment randomly. At low temperatures, the atoms cluster together, i.e. they solidfy. Notice that they naturally for a hexogonally close packed structure (each atom wants to have 6 neighbors). This is not coded into the model anywhere. This is simply the lowest energy configuration do the the shape the Lennard-Jones potential and atoms tend towards it naturally as they try to minimize their energy. 


## THINGS TO TRY

1. Try solving for the minimum energy distance in terms of σ. (if you need a hint on how to do this see the **HINTS** section at the bottom).
2. Try starting the model in the two different configuration options. Then try moving the temperature up and down. See what happens. If you don't start with the HCP structure, are you ever able to get it by adjusting the temperature? Why or why not?
3. Figure out what the lowest energy per particle state is. Why is it this number?
4. See if you can estimate the temparture ranges at which the Lennard-Jones fluid is a solid, liquid and gas (liquid will be difficult to estimate in 2-D fluid).
5. Trying running the model with very few atoms (on the order of 6-10) and then with a lot of atoms (several hundred). Does the temperaturea at which the atoms become a gas (when they start in HCP) change depending on how many there? Why or why not?

## EXTENDING THE MODEL

1. Add a slider to vary EPS (ϵ). See how the temperature ranges at which the atoms are stable as a solid changes as you vary EPS. 
2. Try adding a second type of atom that has a different potential well depth (a different ϵ) when interacting with other atoms of its same type and with atoms of the different type. (ϵ1 could be for type 1 atoms interacting with themselves, ϵ2 could be for type 2 atoms interacting with themselves, and ϵ12 could be for atoms of type 1 and 2 interacting with each other). 
3. Try adding sliders for ϵ1, ϵ2 and ϵ12. See what different types of behavior you can get by varying the potential well depths. 
4. Try making the different atom types be different sizes as well (different σ). There will need to be σ1, σ2 and σ12 (which is (σ1+σ2)/2). 


## RELATED MODELS

Lennard-Jones Chemical Potential
Grand-Canonical 

## HINTS

Hint on solving for the minimum energy: take the derivative of the Lennard-Jones potential and set it equal to zero. 

## CREDITS AND REFERENCES

Original implementation: Jacob Kelter with Prof. Eric Luijten, Northwestern Material Science and Engineering. 

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

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

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

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment1" repetitions="1" runMetricsEveryStep="false">
    <go>go</go>
    <final>file-close-all</final>
    <timeLimit steps="1000000"/>
    <enumeratedValueSet variable="N">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial-Config">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Rho">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 1.0 0.0
0.0 1 1.0 0.0
0.2 0 1.0 0.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
