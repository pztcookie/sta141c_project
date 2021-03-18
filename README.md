# blblm

<!-- badges: start -->
<!-- badges: end -->

## Examples

``` r
library(blblm)

model <- blblm(mpg ~ wt + hp, data = mtcars, m = 3, B = 100)

confint(model, c("wt", "hp"))


sigma(fit, confidence = T)

predict(fit, data.frame(wt = c(3.5, 4), hp = c(160, 180)))

```
