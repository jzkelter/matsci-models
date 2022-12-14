globals [
  electron-diameter
  max-move-dist
  cutoff-dist
  potential-offset
  k0
  pxcor-list
  running-ave-e-distribution
  running-ave-h-distribution
  sorted-meters
  charge-flow
]

breed [electrons electron]
breed [holes hole]
breed [meters meter]

turtles-own [charge]
patches-own [pcharge]
meters-own [
  x-charge
  running-ave-charge
  e-potential
  running-ave-potential
  e-field
  running-ave-field
]

to setup
  clear-all

  ; SET CONSTANTS
  set k0 1
  set potential-offset 0.2
  set pxcor-list sort [pxcor] of patches with [pycor = 0]


  ; CREATE ELECTRONS/HOLES
  set electron-diameter .35
  set cutoff-dist 10 * electron-diameter
  set max-move-dist electron-diameter
  ask patches [
    if-else pxcor < (min-pxcor + world-width / 2)
    [
      set pcolor blue
      if random-float 1 < p-doping% [
        set pcharge -1
        sprout-holes  1[
          hole-init
        ]
      ]
    ]
    [
      set pcolor green - 2

      if random-float 1 < n-doping% [
        set pcharge 1
        sprout-electrons 1 [
          electron-init
        ]
      ]
    ]
  ]

  ; CREATE CONTACTCS
  ask left-contacts [set pcolor grey]
  ask right-contacts [set pcolor grey]

  create-the-meters

  ; SET INITIAL ELECTRON AND HOLE DISTRIBUTIONS
  set running-ave-e-distribution xcor-distribution electrons
  set running-ave-h-distribution xcor-distribution holes

  set charge-flow 0

  reset-ticks
end


to create-the-meters
  ; CREATE METERS
  ask patches with [pycor = 0 ] [
    sprout-meters 1 [
      set hidden? true
    ]
  ]
  create-meters 1 [
    set xcor 0.5
    set hidden? true
  ]
  set sorted-meters sort-on [xcor] meters

end

to electron-init
    set color orange
    set shape "circle"
    set size electron-diameter
    set charge -1
end

to hole-init
  set color black
  set shape "circle"
  set size electron-diameter
  set charge 1
end


;***********************************
;**********GO PROCEDURES************
;***********************************

to go

  generate-electron-hole-pairs

  ask charge-carriers [
    move
    annihilate
  ]
  hold-contacts-charge-neutral
  calc-field

  tick
end

to hold-contacts-charge-neutral
  let left-charge-balance xcor-charge-sum min-pxcor
  (ifelse
    left-charge-balance < 0 [
      create-holes abs (left-charge-balance) [
        hole-init
        set xcor min-pxcor
        set ycor random-ycor
      ]
       update-charge-flow abs left-charge-balance
    ]
    left-charge-balance > 0 [
      ask n-of left-charge-balance holes with [pxcor = min-pxcor] [die]
      update-charge-flow (- abs left-charge-balance)
    ]

    )

  let right-charge-balance xcor-charge-sum max-pxcor
  (ifelse
    right-charge-balance < 0 [
      ask n-of abs(right-charge-balance) electrons with [pxcor = max-pxcor] [die]

      update-charge-flow (- abs right-charge-balance)
    ]
    right-charge-balance > 0 [
      create-electrons right-charge-balance [
        electron-init
        set xcor max-pxcor
        set ycor random-ycor
      ]
      update-charge-flow abs right-charge-balance
    ]
  )
end

to update-charge-flow [amount]
  if ticks > 100 [
    set charge-flow charge-flow + amount
  ]
end

to generate-electron-hole-pairs
  ask patches [
    if random-float 1 < exp(- EHP-formation-energy / temperature) [
      sprout-electrons 1 [
        electron-init
        fd random-float 1
      ]
      sprout-holes 1 [
        hole-init
        fd random-float 1
      ]
    ]
  ]
end


to annihilate
  ifelse breed = electrons [
    let holes-here holes with [distance myself < electron-diameter]
    if (count holes-here) > 0 [
      ask one-of holes-here [die]
      die
    ]
  ] [
    let electrons-here electrons with [distance myself < electron-diameter]
    if (count electrons-here) > 0 [
      ask one-of electrons-here [die]
      die
    ]
  ]

end

to move
  ; STORE OLD VALUE
  let u-old calc-u
  let old-x xcor
  let old-y ycor

  ;MAKE A RANDOM MOVE
  set heading random-float 360
  fd random-float max-move-dist
  let u-new calc-u

  let delta-u u-new - u-old
  ; if potential energy is increased by the move, accept/reject it with prob proportional to temperature
  if (delta-u > 0) and (random-float 1 > exp( - delta-u / temperature) ) [   ; exp( - delta-u / temperature) is the probability the carrier would have the kinetic energy necessary to move up in potential energy this much
    setxy old-x old-y
  ]

end


to-report calc-u
  ; U = k * q1 * q2 / r for two particles
  let u 0

  ask other charge-carriers in-radius cutoff-dist [
    set u u + (k0 * charge * [charge] of myself) / (distance myself  + potential-offset)
  ]

  ask patches in-radius cutoff-dist [
    set u u + (k0 * pcharge * [charge] of myself) / (distance myself + potential-offset)
  ]

  report u + (charge * voltage-at-xcor xcor)

end

to-report voltage-at-xcor [x-cor]
  report voltage * (x-cor - min-pxcor) / world-width
end


to-report charge-carriers
  report (turtle-set electrons holes)
end

to-report left-contacts
  report patches with [pxcor = min-pxcor]
end

to-report right-contacts
  report patches with [pxcor = max-pxcor]
end

;***************************************
;**********E-FIELD AND E-POTENTIAL******
;***************************************

to-report e-field-force [q xcor1 xcor2]
  report q / ((xcor1 - xcor2) ^ 2) * sign (xcor1 - xcor2)
end

to-report sign [num]
  report ifelse-value num > 0 [1] [-1]
end

to-report xcor-charge-mean [x-cor]
  report xcor-charge-sum x-cor / world-height
end

to-report xcor-charge-sum [x-cor]
  report ((sum [charge] of charge-carriers with [abs(xcor - x-cor) <= .5]) + (sum [pcharge] of patches with [pxcor = x-cor]))
end

to calc-field
  calc-charges
  calc-potentials
  calc-e-fields
end

to calc-charges
  ask meters [
    set x-charge xcor-charge-mean xcor
    set running-ave-charge (weighted-running-ave running-ave-charge x-charge)
  ]
end

to calc-e-fields
  ask meters [
    let my-xcor xcor
    set e-field sum [e-field-force x-charge my-xcor xcor] of other meters
    set running-ave-field (weighted-running-ave running-ave-field e-field)
  ]
end

;to calc-potentials
;  let new-weight .1
;  ask meters [
;    let my-xcor xcor
;    set e-potential sum [potential 1 x-charge my-xcor xcor] of other meters
;    set running-ave-potential (1 - new-weight) * running-ave-potential + new-weight * e-potential
;  ]
;end


to calc-potentials
  let last-val 0
  let last-x (min [xcor] of meters) - 1
  foreach sorted-meters
  [m ->
    ask m [
      set e-potential last-val - (e-field * (xcor - last-x))
      set running-ave-potential (weighted-running-ave running-ave-potential e-potential)
      set last-val e-potential
      set last-x xcor
    ]
  ]
end


to-report weighted-running-ave [old-ave new-val]
  let new-weight .05
  report (1 - new-weight) * old-ave + new-weight * new-val
end



;**********************************
;***********Plotting****************
;**********************************
to plot-x-axis
  foreach (list plot-x-min plot-x-max)
  [x -> plotxy x 0]
end

to plot-y-axis
  foreach (list plot-y-min plot-y-max)
  [y -> plotxy 0.5 y]
end


;*************functions for getting distributions (for electron and hole populations plot)***************
to-report xcor-distribution [agent-set]
  let e-dist n-values world-width [[] -> 0]
  foreach [pxcor - min-pxcor] of agent-set
      [px -> set e-dist replace-item px e-dist ((item px e-dist) + 1)]
  report e-dist
end


to-report update-xcor-distribution-ave [ave-xcor-distribution agent-set]
  let new-distribution (list)
  let new-weight .1
  (foreach ave-xcor-distribution (xcor-distribution agent-set)
    [[old-ave new-val] -> set new-distribution lput (weighted-running-ave old-ave new-val) new-distribution]
   )
 report new-distribution
end


to update-ave-hole-distribution
  set running-ave-h-distribution update-xcor-distribution-ave running-ave-h-distribution holes
end

to update-ave-electron-distribution
  set running-ave-e-distribution update-xcor-distribution-ave running-ave-e-distribution electrons
end
@#$#@#$#@
GRAPHICS-WINDOW
246
512
936
771
-1
-1
22.73333333333334
1
10
1
1
1
0
0
1
1
-14
15
-5
5
1
1
1
ticks
30.0

BUTTON
10
115
76
148
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

SLIDER
8
215
180
248
temperature
temperature
0.1
2
1.0
.1
1
NIL
HORIZONTAL

BUTTON
93
115
156
148
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
1

PLOT
225
15
940
135
Potential
NIL
NIL
-10.0
10.0
-2.0
6.0
false
false
"set-plot-x-range min-pxcor max-pxcor" "clear-plot\nset-plot-x-range min-pxcor max-pxcor"
PENS
"Potential" 1.0 0 -5298144 true "" ";foreach sort-on [xcor] meters \n;[ m -> ask m [plotxy xcor e-potential]]\n"
"ave-potential" 1.0 0 -16777216 true "" "foreach sort-on [xcor] meters \n[ m -> ask m [plotxy xcor running-ave-potential]]"
"x-axis" 1.0 0 -4539718 true "" "plot-x-axis"
"y-axis" 1.0 0 -4539718 true "" "plot-y-axis"

PLOT
227
389
938
509
Hole and Electron Populations
NIL
NIL
0.0
5.0
0.0
10.0
false
false
"set-plot-x-range min-pxcor max-pxcor\nset-plot-y-range 0 world-height \n" "clear-plot\nset-plot-x-range min-pxcor max-pxcor\nset-plot-y-range 0 ceiling (world-height * 1.5)"
PENS
"holes" 1.0 0 -16777216 true "" "update-ave-hole-distribution\n  (foreach (range min-pxcor (max-pxcor + 1)) running-ave-h-distribution\n    [[x y] -> plotxy x y]\n  )\n\n;(foreach (range min-pxcor (max-pxcor + 1)) xcor-distribution holes\n;[[x y] -> plotxy x y])"
"electrons" 1.0 0 -955883 true "" "update-ave-electron-distribution\n  (foreach (range min-pxcor (max-pxcor + 1)) running-ave-e-distribution\n    [[x y] -> plotxy x y]\n  )\n;(foreach (range min-pxcor (max-pxcor + 1)) xcor-distribution electrons\n;[[x y] -> plotxy x y])"
"pen-2" 1.0 0 -4539718 true "" "plot-y-axis"
"holes-instantaneous" 1.0 0 -7500403 true "" ";(foreach (range min-pxcor (max-pxcor + 1)) xcor-distribution holes\n;[[x y] -> plotxy x y]\n;)"
"pen-4" 1.0 0 -2674135 true "" ""

PLOT
226
263
938
383
Charge Density
NIL
NIL
0.0
0.0
-1.0
1.0
true
false
"" "clear-plot\nset-plot-x-range min-pxcor max-pxcor"
PENS
"ave-charge-density" 1.0 0 -16777216 true "" "foreach sort-on [xcor] meters \n[ m -> ask m [plotxy xcor running-ave-charge]]"
"x-axis" 1.0 0 -4539718 true "" "plot-x-axis"
"y-axis" 1.0 0 -4539718 true "" "plot-y-axis"

PLOT
225
138
939
258
E-Field
NIL
NIL
0.0
10.0
-5.0
1.0
false
false
"set-plot-x-range min-pxcor max-pxcor" "clear-plot\nset-plot-x-range min-pxcor max-pxcor\nset-plot-y-range round((min [running-ave-field] of meters) - 1) 3"
PENS
"running-ave" 1.0 0 -16777216 true "" "foreach sort-on [xcor] meters \n[ m -> ask m [plotxy xcor running-ave-field]]"
"x-axis" 1.0 0 -4539718 true "" "plot-x-axis"
"y-axis" 1.0 0 -4539718 true "" "plot-y-axis"

SLIDER
5
18
177
51
p-doping%
p-doping%
0
1
1.0
.1
1
NIL
HORIZONTAL

SLIDER
6
67
178
100
n-doping%
n-doping%
0
1
1.0
.1
1
NIL
HORIZONTAL

SLIDER
8
173
180
206
EHP-formation-energy
EHP-formation-energy
1
15
15.0
1
1
NIL
HORIZONTAL

SLIDER
11
285
183
318
voltage
voltage
-50
50
-240.0
1
1
NIL
HORIZONTAL

PLOT
11
348
211
498
Current
NIL
NIL
0.0
10.0
-1.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "if ticks > 0 [plotxy ticks (charge-flow / (ticks - 99))]"
"x-axis" 1.0 0 -7500403 true "" "plot-x-axis"

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
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
NetLogo 6.1.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1000"/>
    <metric>charge-flow</metric>
    <metric>ticks</metric>
    <enumeratedValueSet variable="n-doping%">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temperature">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-doping%">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="voltage" first="-400" step="20" last="-220"/>
    <enumeratedValueSet variable="EHP-formation-energy">
      <value value="15"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
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
0
@#$#@#$#@
