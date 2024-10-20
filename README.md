# Phillips-Snyder Lifting-Line Theory (PS-LLT)
Phillips, W. F., and Snyder, D. O., "Modern Adaptation of Prandtl's Classic Lifting-Line Theory," Journal of Aircraft, Vol. 37, No. 4, 2000, pp. 662-670.

`DragOpt` does inverse design to create an elliptical span load. The induced downwash (and hence induced angle of attack $\alpha_i$) is constant along the span. From a sectional lift model, a twist distribution can be generated and exported alongside other information about the wing. A portion of code at the end of the program exports the optimized wing in a format that is readily usable with `importWing`. Users can import custom wings using `importWing` by supplying a file that adheres to the following:

| Element | Description |
| --- | --- |
| Delimiter | Any non-numeric character excluding `.` and `\n` |
| Header | (optional) Text |
| Column 1 | x-coordinates of the leading edge defined at "posts" |
| Column 2 | x-coordinates of the trailing edge defined at "posts" |
| Column 3 | Span locations of "posts" |
| Column 4 | z-coordinates of the c/4 line defined at "posts" |
| Column 5 | Twist angles (in degrees) of the sections defined at span locations midway between "posts" with any numeric value padding the bottom of the column |
| Column 6 | Sectional lift-curve slopes (in change per radian) of the sections defined at span locations midway between "posts" with any numeric value padding the bottom of the column |
| Column 7 | Zero-lift angle of attack (in degrees) of the sections defined at span locations midway between "posts" with any numeric value padding the bottom of the column |

The term "posts" refers to nodes in the sense of posts and fences. Columns with values defined at control points (midway between "posts") have one fewer value than columns with values defined at the nodes and therefore need to be padded with a filler value. **The top-to-bottom ordering must be from the wing root to the wing tip. As such, the file must only describe half of a symmetric wing.**

## DragOpt
### Induced Angle of Attack
$$L'=\rho V_\infty \Gamma$$
$$dL=\rho V_\infty \Gamma dy$$
$$L=\sum_{i=1}^N dL_i=\rho V_\infty \sum_{i=1}^N \Gamma_i dy_i=\rho V_\infty b \sum_{i=1}^N \Gamma s_i$$
$$C_L=\dfrac{L}{\frac{1}{2}\rho V_\infty^2 S}=\dfrac{2b\sum_{i=1}^N \Gamma_i s_i}{V_\infty S}=\dfrac{2\sum_{i=1}^N \Gamma_i s_i}{V_\infty c_\text{avg}}$$
$$C_L=\sum_{i=1}^N g_i s_i$$
$$g=\dfrac{2\Gamma}{V_\infty c_\text{avg}}$$
$$\alpha_i=\dfrac{\Gamma_0}{2bV_\infty}=\dfrac{g_0 c_\text{avg}}{4b}$$