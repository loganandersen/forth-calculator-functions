: squared ( n -- n^2 )
    dup * ;

: .p ( a -- a )
    dup . ;

: f.p ( r -- r )
    fdup f. ;

: fsquared ( r -- r^2 )
    fdup f* ;

: fcubed ( r -- r^2 )
    fdup fsquared f* ; 

\ euclidean distance 2d
: fdist2 ( r1 r2 -- r3 )
    fsquared fswap
    fsquared f+ fsqrt ;

: frev3 ( r1 r2 r3 -- r3 r2 r1 )
    fswap frot ;

\ : frev4 ( r1 r2 r3 r4 -- r4 r3 r2 r1 )

\ first argument is the x value, then coefficents a b c
: fquadratic ( rx ra rb rc -- rs )
    \ reverse stack a in front 
    frev3
    \ a x^2
    3 fpick fsquared f*
    \ b x 
    fswap 3 fpick f*
    \ ax^2 + bx + c and drop x
    f+ f+ fnip ;

: fradtodeg ( r -- rd )
    360e pi 2e f* f/ f* ;

: fdegtorad ( rd -- r )
    2e pi f* 360e f/ f* ;

: fdegcos ( rtd -- r )
    fdegtorad fcos ;

: fdegsin ( rtd -- r )
    fdegtorad fsin ;

: fdegtan ( rtd -- r )
    fdegtorad ftan ;


: fa,b*c ( ra rb rc -- rac rbc )
    fdup frot ( a c c b )
    f* frot frot ( bc c a )
    f* fswap ( ac bc ) ; 

: fatan2deg ( rsin rcos -- rtheta )
    fatan2 fradtodeg ;

: fxytodeg ( rx ry -- rt )
    fswap fatan2deg ;

: fxydegmagnitude ( rm rt -- ry rx )
    fdegtorad fsincos ( mag sin cos )
    frot ( sin cos mag )
    fa,b*c ;

\ moments of intertia table
\ Slender rod (center)          : I = 1/12 M L^2
\ Hollow Cylender               : I = 1/2 M(R_1^2 + R_2^2)
\ Solid cylender                : I = 1/2 MR^2
\ ring cylinder                 : I = MR^2 \these can be found from the hollow cylendar formula
\ Solid Sphere                  : I = 2/5 MR^2
\ Hollow Sphere                 : I = 2/3 MR^2
\ Rectangular plate center      : I = 1/12 M R^2
\ thin Rectangular plate edge   : I = 1/3 M R^2

: solid-cylinder-moment ( mass radius -- moment-of-inertia )
    fsquared 0.5e f* f* ;

: square-2-args ( a b -- a^2 b^2 ) 
    fsquared fswap fsquared fswap ;

: square-2-and-add ( a b -- a^2+b^2 )
    square-2-args f+ ;

: plate-moment ( mass length-a length-b -- moment-of-intertia )
    square-2-and-add 1e 12e f/ f* f* ;

: hollow-cylinder-moment ( mass rad1 rad2 -- moment of intertia )
    square-2-and-add 0.5e f* f* ;
    
: cm->m ( cm -- m )
    \ centemeters to meters
    100e f/ ;

: slender-rod-moment ( mass length -- moment-of-inertia )
    fsquared f* 1e 12e f/ f* ;