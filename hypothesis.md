Hypothesis (probably false): You can use a "bag-of-words" model to compute the probability that a bound should break, and low probability bonds will be broken in symmetrical molecules.

  - HOWTO: Compute the probability of each word for a bond; then, take the product.  This is not numerically stable, so try:
  	
  		Sum(- log [word probability]) = - log (Product [word probability])
 
  	This is the *naïve Bayesian* likelihood.
  	
  - Where to get words:
     - Matt's approach: Take 200 molecules, and pair them up randomly.  This gives you (200 choose 2) = 19900 combinations. 	*5000 choose 2 = 12497500
      Then, compute the *Maximum Common Substructure* for these 19900 and uniquify; this gives you ~3000 unique common substructures
        - Try on common molecules: Matt predicts even fewer MCSs
        - Does not preserve chirality
        
     - Prof. G's approach: Try natural products
        - Hypothesis: more = better (Matt disagrees)
   
   - How to compute frequencies:
      - Use `hasSubstructureMatch` on some number of molecules
      
Supervised Learning Approach:
   - Observation: any algorithm that suggests what bonds to cut based on a bag-of-words model will suggest symmetrical cuts with equal weight
   - 