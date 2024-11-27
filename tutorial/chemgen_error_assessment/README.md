# Chemgen


## Table of Contents

- [Generating Custom Test](#Generating Custom Test)
- [Generating Error Data](#Generating Data)
- [Post Processing Data](#Post Processing)

## Generate Custom Test
In all tutorials we will shorten the use of ChemGen's paths to the repo. To access ChemGen in this folder run the following command:

```bash
export PATH="$(cd ../../bin && pwd):$PATH"
```

Now, ChemGen can be run from any directory by simply specifying `chemgen.py`

ChemGen has an option "--custom-test" where a python function, `write_test` can be overwritten to create a custom `chemgen.cpp`. We've included `custom_test.py` in this directory.


## Generating Data

#### \(10^{-4}\) (Excellent Accuracy)
- **Meaning**: The approximate value is extremely close to the true value.
- A log-scale error of \(10^{-4}\) corresponds to a deviation factor:
  \[
  10^{10^{-4}} \approx 1.0002303
  \]
- This is essentially a **tiny fractional difference**, less than **0.023%**.

---

#### \(10^{-1}\) (Moderate Accuracy)
- **Meaning**: The approximate value is moderately close to the true value.
- A log-scale error of \(10^{-1}\) corresponds to a deviation factor:
  \[
  10^{10^{-1}} \approx 1.2589
  \]
- This means the approximate value is within about **25.89%** of the true value. 
- Such deviations are often acceptable in systems with high variability or where perfect precision is not required.

---

### Summary of Accuracy Levels

| Log-Scale Error | Deviation Factor \(10^{\text{Log-Error}}\) | Description                 |
|------------------|-------------------------------------------|-----------------------------|
| \(10^{-4}\)      | \(1.0002303\)                            | Excellent accuracy (<0.023%) |
| \(10^{-1}\)      | \(1.2589\)                               | Moderate accuracy (~25.89%) |

---

This formatting ensures mathematical clarity and makes the content easy to read and understand in Markdown-rendering tools.

## Post Processing
