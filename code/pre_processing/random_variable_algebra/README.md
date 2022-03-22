# Random variable algebra with `pacal`

Generating pdf and cdf of the product of random variables
  - reliability is governed by the tails, so we have to accurately estimate those
  (simulation based approaches are inefficient)
  - pacal implement more efficient algorithms which allow for efficiently performing
    algebra of random variables
  - the outputs serve as input in the reliability analysis

## Running the code

Install dependencies:

```commandline
pip install -r requirements.txt
```

Run `pacal_product.py`. Results are saved in the `results` folder.

## Development notes

* The need for the product distribution is due to the design decisions made in the 
calibration code (back then it seemed to be easier to create these products.
rather than restructuring the code though this decision could and should be revisited)
* Developed under Windows 10, for the used package version numbers see 
`requirements_dev.txt`