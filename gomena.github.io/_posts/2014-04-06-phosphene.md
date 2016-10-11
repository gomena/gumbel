---
layout: post
categories: gifs
author: Gonzalo Mena
excerpt_separator: <!--more-->
comments: true
---

Have you ever seen shapes or patterns when you close your eyes and rub your eyes? This very personal phenomenon of experiencing vision with no light is called 'phosphene', and is thought to be caused by non-photonic stimulation of neurons in the retina. This GIF is just an attempt to communicate this experience. How I did it? Basically it is a wave that propagates from the center. At each frame colors are updated in regions where the wave front is active, and the rule to update the color is to just add a constant (mod 1) to the color variable, which is later translated into an actual color using a predefined colormap.

![](http://stat.columbia.edu/~gonzalo/gallery/CirclePixels.gif)

# Mandala
Mandalas are pictorial representations of the universe, according to hinduist and buddhist cosmogonies. When doing this other GIF I realized that when changing the color update rule by the rule given by the ergodic transformation T(x)=x+π (the key is that π is an irrational number) I obtained this other figure. So you can think of this as a moving Mandala by using the above transformation. Note: Maybe there is no sound theoretical basis to justify the use of that transformation, but at least gives me a reason to use this appealing title


![](http://stat.columbia.edu/~gonzalo/gallery/CirclePicompressed.gif)
