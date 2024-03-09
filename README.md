# Predator Prey Pursuits

**Context:**  
During my master's in integrative biology, I took a course called ecological and behavioral modeling, where we learned how to implement various biological models in R. This course culminated with a self-guided project. In this project, I modified an existing model to find the optimal prey flight strategy based on a predator's pursuit strategy. The following sections outline the derivation of the model, its implementation, and my findings.

**Introduction:**  
Predator-prey pursuits dwell at the intersection of biomechanics and evolutionary biology. The outcomes of these pursuits determine the survival of both organisms involved - the predator by acquiring or failing to acquire resources and the prey by being captured or not. Thus, adaptations that aid or harm these pursuits will be selected for or against, respectively (Irschick et al. 2008). Understanding more about predator-prey pursuits can give insight into the dynamics of the selective pressures that act on organisms.
Predator pursuit strategies can vary widely, ranging from endurance-exhaustion to high-speed short-term (Guinet et al. 2007). It is likely that the prey item’s strategy will depend on the strategy that its respective predator takes. Previous research has focused primarily on the predator’s perspective (Wilson et al. 2013, 2015). The papers that do investigate the prey’s perspective focused on turning dynamics, were taxonomically narrow, or ignored optimizing prey strategy (Weihs and Webb 1984, Hammond et al. 2020).
To alleviate these issues, I used a modified version of a previously developed model to investigate the optimal strategy of a generalized prey item. This prey could choose between different strategies. I hypothesize that the best prey strategy choice will be to match the predator’s strategy.

**Predator Model:**  
In their 1977 paper, Elliot et al. used the following model for predator-prey pursuits in various African taxa (Elliott et al. 1977): 
  
$\frac{dx}{dt}=a(1-e^{-bt})$									(1)  
  
The equation models one-dimensional velocity as a function of time, where the parameters $a$ and $b$ control the maximum velocity and acceleration, respectively. This equation produces velocity-time curves that asymptote to the maximum velocity (Figure 1). Infinitely sustaining a velocity, though theoretically possible for organisms within their aerobic capacity, is unrealistic in most cases (Bassett and Howley 1997). To alleviate this, I added linear decay to the model, producing the following equation:  
  
$\frac{dx}{dt}=a(1-e^{-bt})-cx$								(2)	  
  
(2) acts similarly to (1), but the addition of $cx$ controls the rate of fatigue. Higher values of $c$ produce quicker fatigue but also significantly decrease the maximum velocity. Given that I wanted to modify each parameter independently, I modified (2) to keep the maximum velocity constant $∀ c$. The linear decay causes a discrete maximum velocity rather than an asymptote (Figure 2). Using the first derivative test, I derived the following value for the time at which this occurs:  
  
$t_{max}=\ln⁡(\frac{ab}{c}) b^{-1}$   									(3)  
  
I then multiplied (2) by an unknown constant $a^*$ such that $\frac{dx}{dt}$ equals $a$ at $t_{max}$. Substituting these into (2) produces the following:  
  
$a= a^* (a(1-e^{\frac{-b}{b}  ln⁡(\frac{ab}{c}) })-cx)$  
    $= a^* (a-\frac{a^2 b}{c-cx})$									(4)  
      
Solving (4) for $a^*$ and substituting produces the following:  
  
$\Large\frac{dx}{dt}=\frac{(a^2 (1 - e^{-bt} )  - act)}{(a - \frac{c}{b(1 + ln⁡(\frac{ab}{c}))})}$								(5)  
  
(5) models the velocity of an animal, where $a$ controls the animal’s maximum velocity, $b$ controls its acceleration – higher values produce faster acceleration – and $c$ controls its rate of fatigue after reaching its max velocity – higher values produce quicker fatigue. It is important to note that $\frac{ab}{c} ≠ 1 \text{ and } \frac{ab}{c} > 0$ due to the logarithm in the denominator. This restricts the values that can be tested.  
  
I then substituted $t$ for a variable $t^*$ and integrated from 0 to an arbitrary time $t$ to produce the following:  
  
$\large x(t) \Large =  \frac{ a^2 t + a^2 b(e^{-bt}-1)  - \frac{a}{2} ct^2}{a - \frac{c}{b(1 + ln⁡(\frac{ab}{c}))}}$							(6)  
   
(6) models the position of an animal, moving according to (5), at any given time t. I used this to model the movement of the predator. 
   
 ![image](https://github.com/maldridge3/Predator_Prey_Pursuits/assets/96790192/6376ecdb-c0a1-41b1-b824-fe5805490300)  
**Figure 1:** This shows the velocity vs time curves produced by the original model. Note that time is in seconds, in contrast with the model used in the current paper. (Figure from Elliot et al. 1977)

![image](https://github.com/maldridge3/Predator_Prey_Pursuits/assets/96790192/15bd2894-c1b8-4920-bdf3-c23b4931c0fd)  
**Figure 2:** This shows the velocity vs time curves for the three locomotive strategies. Green is type 1, red is type 2, and orange is type 3. The shape of the curves results in a single maximum. Time is not necessarily in seconds. Note that type 3 also fatigues, but this is not noticeable until larger times. (Figure produced using Geogebra)

**Prey Model:**  
I made several modifications to (6) for the prey’s model. First, I added a location parameter $r$ to control the start time. This is the prey’s reaction time, with $r > 0$ corresponding to a predator-initiated pursuit and $r < 0$ corresponding to a prey-initiated pursuit. I also added a parameter $d$ to control the prey’s starting distance from the predator. Adding these to (6) produces the following:  
  
$\large x(t) = d + \Large \frac{a^2 (t-r) + a^2 b(e^{-b(t-r)}-1)  - \frac{a}{2} c(t-r)^2}{a - \frac{c}{b(1 + ln⁡(\frac{ab}{c}))}}$				(7)  
  
Lastly, I added a parameter $M$ (mod coefficient) to modify the graph about the base strategies. To derive this, I let the prey’s acceleration $b_{prey}$ be a multiple of the predator’s acceleration such that $b_{prey} = M·b_{pred} = M·b$. I did the same for the fatigue coefficient $c_{prey}$ such that $c_{prey} = M·c_{pred} = M·c$. I did not do this for the maximum velocity, as the goal was to alter only the shape of the velocity curve about a baseline. Substituting $b_{prey}$ and $c_{prey}$ into (7) produces the following:  
  
$\large x(t)=d+ \Large \frac{a^2 (t-r) + a^2 Mb(e^{-Mb(t-r)}-1)  - \frac{a}{2} Mc(t-r)^2}{a - \frac{c}{b(1 + ln⁡(\frac{ab}{c}))}}$			(8)  
  
(8) models the position of the prey item at any time $t$. Note that a mod coefficient $M = 1$ does not alter the function from the baseline strategy, mod coefficients $M ∈ (0,1)$ shift the prey’s velocity curve towards high speed and acceleration, decreasing endurance, and mod coefficients $M > 1$ shift the prey’s velocity towards endurance, decreasing acceleration and max velocity. 

**Assumptions:**  
Since the model is one-dimensional, it assumes that the prey flees directly away from the predator. This applies in several scenarios. First, since larger organisms cannot turn as quickly as slower organisms, prey that are larger than their predators will generally flee in a straight line (Wilson et al. 2015). Second, at larger distances, prey will flee in a straight line to prevent predators from cutting corners (Wilson et al. 2015). Lastly, in some predator-initiated pursuits, prey will flee in a straight line (Nishiumi and Mori 2020). 
I assume that most predator-prey locomotive strategies can be categorized into three types. Type 1 is an endurance-exhaustion strategy. Type 3 is high speed and acceleration at the cost of increased fatigue. Lastly, type 2 is an intermediate strategy. These can be viewed as an organism’s choice of velocity with respect to its metabolic capacity. Higher velocities and accelerations lead to less sustainable bouts.
The model also assumes that the predator and prey have similar aerobic capacities. That is to say: what is considered a sprint for the predator is also considered a sprint for the prey. I assume that this relationship is all that matters. The absolute velocities of each strategy are not important. Because of this, the times and distances produced by this model cannot be related to physical measurements.
I assumed the prey escaped in three situations. First, if the prey exceeded a certain radius, I assumed the predator would either lose the prey or give up. Second, the prey escaped if the model reached 180 timesteps. This primarily applied to type 3 predators, as the other two never reached this timestep. Third, the prey escaped if the predators reached 50% of their max velocity. This was to better reflect the reality that organisms generally do not run until they physiologically cannot move anymore. In contrast to the second criterion, the third primarily applied to type 1 and 2 predators, as time would run out long before a type 3 predator reached this velocity.
I assumed the prey was caught in two situations. First, if the distance between the predator and the prey was less than 0.1, the predator caught the prey. The radius was necessary due to the models being evaluated at discrete timesteps (e.g. their positions were never exactly equal). Second, the prey was captured if the predator passed the prey between timesteps. This was to avoid scenarios where the predator was much faster than the prey and passed the prey between timesteps. I assumed that a predator this fast could catch the prey.

**Methods:**  
I utilized the R programming environment to run my model. The program tested each predator prey combination for mod coefficients $M ∈ [0.6,1.4]$ in intervals of 0.2, reaction times $r ∈$ $[-0.5,0.5]$ in intervals of 0.25, and starting distances $d ∈ [1,10]$ in intervals of 1. The program progressed in timesteps of 0.1 (note that this is not necessarily in seconds). R plugged the time into each function and compared the predator and prey positions. If the positions satisfied the capture or escape conditions, the program logged the result with its corresponding conditions into a matrix.   
  
Next, the program ordered the prey types in order of success for each predator type. Then, it computed the marginal survival rates for the mod coefficients, reaction times, and starting distances for first and second most successful prey types (note that using marginal rates sometimes produces survival rates less than 1, even though the model is deterministic). These are called “optimal” and “sub-optimal”, respectively. The program finished by producing graphs of these data. I repeated this several times for various radii of escape and relative velocities for type 1 and 3 strategies to see how these influenced the results.

**Results:**  
Matching predator strategy was the optimal strategy for predator type. Type 2 was sub-optimal for type 1 and 3 predators, while type 1 was sub-optimal for a type 2 predator. It is important to note that all three prey types were viable against a type 3 predator. Increasing the difference in velocity between type 1 and type 3 decreased the survival of type 1 and 2 prey against a type 3 predator. 
The radius of escape only slightly influenced prey survival against a type 3 predator. Type 1 and 2 prey remained viable against a type 1 predator for moderate radii of escape, though only type 1 prey were viable at high radii of escape. Type 1 and 3 prey were only viable against a type 2 predator at low radii of escape. 
    
![image](https://github.com/maldridge3/Predator_Prey_Pursuits/assets/96790192/bab6cfdd-dd3b-4307-b0f8-20e38526e7cd)  
**Figure 3:** This shows prey survival (y-axis) for each predator type (1-3, top to bottom) against mod coefficient (x-axis). The color shows which prey type was optimal for each predator (type 1 = peach, type 2 = green, type 3 = blue).
  
![image](https://github.com/maldridge3/Predator_Prey_Pursuits/assets/96790192/d0df9470-3550-4dc1-bec7-9712945f62d9)  
**Figure 4:** This shows optimal prey survival (y-axis) against reaction time (x-axis) for each predator type (1-3, top to bottom). The color shows which prey type was optimal for each predator (type 1 = peach, type 2 = green, type 3 = blue).

Mod coefficients did not affect survival in optimal prey types (Figure 3). This was also true for sub-optimal prey (not shown due to redundancy). Survival decreased with reaction time (Figure 4). Interestingly, sub-optimal prey could still survive at high rates so long as the pursuit was prey-initiated (Figure 5). Prey survival increased with starting distance (Figure 6). This was correlated more strongly in type 3 predators (type 1: ρ ≈ 0.78; type 3: ρ ≈ 0.94). 

![image](https://github.com/maldridge3/Predator_Prey_Pursuits/assets/96790192/94486bb2-0463-4533-8a75-d0f8a9c45a23)  
**Figure 5:** This shows sub-optimal prey survival (y-axis) against reaction time (x-axis) for each predator type (1-3, top to bottom). The color shows which prey type was sub-optimal for each predator (type 1 = peach, type 2 = green, type 3 = blue).
  
![image](https://github.com/maldridge3/Predator_Prey_Pursuits/assets/96790192/d4e2b820-52e2-4869-a61b-f72dacade628)  
**Figure 6:** This shows optimal prey survival (y-axis) against starting distance (x-axis) for each predator type (1-3, top to bottom). The color shows which prey type was optimal for each predator (type 1 = peach, type 2 = green, type 3 = blue).


**Discussion:**  
The results support my hypothesis that prey should match their predator’s strategy. However, there are two caveats: 
First, though type 3 prey were most successful against a type 3 predator, there was only about a ten percentage-point difference between the most and least successful prey types. It is possible that this is simply because type 3’s velocity was not high enough in relation to other strategies, as the difference between optimal and sub-optimal survival rates increased with the difference between type 1 and 3 velocities. In addition, the areas under the velocity curves are similar at early time steps with the chosen relative velocities. This implies that prey were covering similar distances in the early timesteps. Follow-up research investigating what a reasonable relative velocity choice should be is needed to determine if this result is an artifact of experimental design. 
Second, though sub-optimal prey types were generally not as successful, they were reasonably successful so long as the pursuits were prey-initiated. Thus, should there be a selective pressure that requires adaptations causing sub-optimal locomotive performance, prey can adapt to this pressure and remain successful in predator-prey pursuits so long as there is a concomitant adaptation to their senses. This allows them to spot predators at larger distances and makes them more likely to initiate the pursuit before the predator does, both of which increased survival regardless of prey strategy.
It is interesting that starting distance is correlated more strongly in type 3 predators. This is likely due to the greater importance of the initial pursuit. A type 3 predator is going to be faster than most of the prey strategies pitted against it. Thus, extra distance can help slower prey stay outside the radius of capture until the predator fatigues. In contrast, a type 3 predator is going to be slower than most of the prey it encounters. Thus, extra distance is only helpful if it helps the prey get outside the escape radius.
Perhaps the most interesting result was the lack of any effect by the mod coefficients. This implies that there is flexibility about any given strategy. Should adaptations with minor effects on prey locomotion be desirable, prey can realize these without negatively impacting success in predator-prey pursuits. However, it should be noted that adding time-dependent capture success may affect this result. Previous studies modeled the probability of capture success decreasing as pursuit time increased (Elliott et al. 1977). Further analysis of my own results showed that, though the outcome of pursuits did not change with different mod coefficients, the distances covered did change. Thus, the prey may experience higher escape rates at different mod coefficients under a time-dependent probabilistic model. Further research investigating how mod coefficients affect survival in probabilistic conditions is needed.


**Conclusions:**  
Here, I presented a modification of a previous model that allows greater control of the shape of predator-prey velocity curves. I used this model to investigate the dynamics of predator-prey pursuits. I defined three locomotive strategies based on the acceleration, max velocity, and endurance of an organism. I hypothesized that prey should match the locomotive strategy of their predators and found support for this hypothesis. Though my results suggest there is flexibility around ideal strategies, future research is needed to verify this. This research provides insight into the evolution of predator-prey pursuits and illuminates how organisms can adapt to selective pressures.

**References:**  
Bassett, D. R. J., and E. T. Howley. 1997. Maximal oxygen uptake: “classical” versus “contemporary” viewpoints. Medicine & Science in Sports & Exercise 29:591–603.  
Elliott, J. P., I. M. Cowan, and C. S. Holling. 1977. Prey capture by the African lion. Canadian Journal of Zoology 55:1811–1828.  
Guinet, C., P. Domenici, R. de Stephanis, L. Barrett-Lennard, J. K. B. Ford, and P. Verborgh. 2007. Killer whale predation on bluefin tuna: exploring the hypothesis of the endurance-exhaustion technique. Marine Ecology Progress Series 347:111–119.  
Hammond, J. E., S. Witkowski, T. Wilson, C. A. Zouvi, N. L. Goetz, N. F. Eck, and R. W. Clark. 2020. Know Thine Enemy: Predator Identity Influences the Response of Western Banded Geckos (Coleonyx variegatus) to Chemosensory Cues. Journal of Herpetology 54:480–484.  
Irschick, D. J., * Jerry J. Meyers, J. F. Husak, and J.-F. L. Galliard. 2008. How does selection operate on whole-organism functional performance capacities? A review and synthesis. Evolutionary Ecology Research 10:177–196.  
Nishiumi, N., and A. Mori. 2020. A game of patience between predator and prey: waiting for opponent’s action determines successful capture or escape. Canadian Journal of Zoology 98:351–357.  
Weihs, D., and P. W. Webb. 1984. Optimal avoidance and evasion tactics in predator-prey interactions. Journal of Theoretical Biology 106:189–206.  
Wilson, A. M., J. C. Lowe, K. Roskilly, P. E. Hudson, K. A. Golabek, and J. W. McNutt. 2013. Locomotion dynamics of hunting in wild cheetahs. Nature 498:185–189.  
Wilson, R. P., I. W. Griffiths, M. G. Mills, C. Carbone, J. W. Wilson, and D. M. Scantlebury. 2015. Mass enhances speed but diminishes turn capacity in terrestrial pursuit predators. eLife 4:e06487.





