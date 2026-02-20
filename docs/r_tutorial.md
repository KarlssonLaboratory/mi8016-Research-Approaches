# Introduction to R: A Beginner's Guide

## Table of Contents

1. [Who This Guide Is For](#who-this-guide-is-for)
2. [Installing R](#installing-r)
3. [Your First R Commands](#your-first-r-commands)
4. [Choosing an IDE](#choosing-an-ide)
   - [RStudio](#rstudio)
   - [Positron](#positron)
5. [R Packages](#r-packages)
6. [Mini Project: Your First Data Analysis](#mini-project-your-first-data-analysis)
7. [Troubleshooting](#troubleshooting)
8. [Getting Help](#getting-help)

---

## Who This Guide Is For

This tutorial is for anyone who is new to R and wants to get up and running quickly — no prior programming experience required. By the end, you will have R installed, know how to run basic commands, and have completed a simple data analysis project using real data.

---

## Installing R

[R](https://www.r-project.org/) is a free ([as in freedom!](https://www.gnu.org/philosophy/free-sw.html)), open-source environment for statistical computing and graphics. It runs on Windows, macOS, and Linux.

The R environment is a fully planned and coherent system, rather than an incremental accretion of very specific and inflexible tools, as is frequently the case with other data analysis software (e.g. GraphPad Prism, IBM SPSS).

Download R from the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/), which hosts up-to-date versions of R and its documentation:

- [R for Windows](https://cran.r-project.org/bin/windows/base/)
- [R for macOS](https://cran.r-project.org/bin/macosx/)
- [R for Linux](https://cran.r-project.org/bin/linux/)

Once installed, you can launch R from the command line (Linux/macOS) or open the R application on Windows. However, most users work with an Integrated Development Environment (IDE) rather than the base interface — more on that in the [Choosing an IDE](#choosing-an-ide) section.

---

## Your First R Commands

Once R is installed, open it and try the following commands to get comfortable with the basics.

### Basic Arithmetic

R works like a calculator. Type a command after the `>` prompt and press `Enter`:

```r
> 2 + 2
[1] 4

> 10 * 5
[1] 50

> 100 / 4
[1] 25
```

The `[1]` before the result is R's way of labeling the first element of the output — you can ignore it for now.

### Storing Values in Objects

You can save values using the assignment operator `<-` and reuse them later:

```r
a <- 2 + 2    # stores 4
b <- 5 * 2    # stores 10
c <- 12 / 3   # stores 4
```

Call the variable name to print its value:

```r
> a
[1] 4
> a + b + c
[1] 18
> a - b
[1] -6
```

### Working with Vectors

A vector is a sequence of values — one of R's most fundamental data structures:

```r
# Create a vector of numbers
ages <- c(23, 35, 42, 28, 19)

# Compute summary statistics
mean(ages)    # average
median(ages)  # middle value
sum(ages)     # total
```

### Saving and Running Scripts

Instead of typing commands one by one, you can save them in a script — a plain text file with the `.R` extension. To run a script from the command line:

```bash
Rscript my_script.R
```

In an IDE, you can run a script line by line with `Ctrl + Enter`, or all at once using the **Source** button.

---

## Choosing an Integrated Development Environment (IDE)

While R comes with a basic graphical interface, most users find a dedicated IDE much more productive. Here are the two most popular options.

### RStudio

[RStudio](https://posit.co/download/rstudio-desktop/) is the most widely used IDE for R and the usual choice for beginners. Its interface is organized into four panes:

1. **Source** — write and save your R scripts here
2. **Console** — run commands interactively
3. **Environment** — see the variables and data you've created
4. **Output** — view plots, tables, files, and help documentation

![RStudio Panes](images/rstudio-panes-labeled.jpeg)

You can run code line by line with `Ctrl + Enter`, or run the entire script using the **Source** button in the top-right of the Source pane. Detailed guides are available at [docs.posit.co](https://docs.posit.co/ide/user/ide/get-started/).

### Positron

[Positron](https://positron.posit.co/download.html) is a newer IDE from the same team that built RStudio. It is built on the open-source foundation of Visual Studio Code and is a good fit if you:

- work with multiple languages (Python, R, Bash, Quarto, etc.)
- want access to the VS Code extensions marketplace
- prefer a more customizable environment

For beginners focused on R, RStudio is perfectly sufficient. Positron becomes more appealing as your work expands to other languages or tools.
However, some data analysis features of Positron makes it even better than RStudio.

![RStudio vs Positron comparison](images/user-interface-rstudio-vs-positron.jpeg)

---

## R Packages

R's functionality can be extended with **packages** — collections of functions, data, and documentation contributed by the community. There are thousands available.

### Installing from CRAN

[CRAN](https://cran.r-project.org/web/packages/available_packages_by_name.html) is the main repository for R packages. Install a package with:

```r
install.packages("ggplot2")  # a popular data visualization package
```

You only need to install a package once. To use it in a session, load it with `library()`:

```r
library(ggplot2)
```

### Installing from Bioconductor

Packages for bioinformatics are typically hosted on [Bioconductor](https://bioconductor.org/packages/release/BiocViews.html#___Software) rather than CRAN. To install them, first install the `BiocManager` package:

```r
install.packages("BiocManager")
```

Then use it to install Bioconductor packages:

```r
BiocManager::install("GenomicFeatures")
```

The `::` notation specifies which package a function comes from — useful when two loaded packages share a function name.

---

## Mini Project: Your First Data Analysis

Let's put everything together with a short end-to-end analysis using the built-in `iris` dataset. Collected by botanist Edgar Anderson in 1935 and made famous by statistician Ronald Fisher, it contains measurements (in centimeters) of sepal length, sepal width, petal length, and petal width for 150 flowers across three species of the genus Iris: *setosa*, *versicolor*, and *virginica*. It's one of the most widely used datasets in data science and comes built into R — no downloading required.

![Iris species](images/iris.png)

### Step 1: Explore the Data

A good way of starting to explore the data is using the `str()` function from the `utils` package that comes with base R and it is automatically loaded when you start a session.

```r
str(iris)
'data.frame':	150 obs. of  5 variables:
 $ Sepal.Length: num  5.1 4.9 4.7 4.6 5 5.4 4.6 5 4.4 4.9 ...
 $ Sepal.Width : num  3.5 3 3.2 3.1 3.6 3.9 3.4 3.4 2.9 3.1 ...
 $ Petal.Length: num  1.4 1.4 1.3 1.5 1.4 1.7 1.4 1.5 1.4 1.5 ...
 $ Petal.Width : num  0.2 0.2 0.2 0.2 0.2 0.4 0.3 0.2 0.2 0.1 ...
 $ Species     : Factor w/ 3 levels "setosa","versicolor",..: 1 1 1 1 1 1 1 1 1 1 ...
```

The `$` operator in R is used to extract a named element from a list or data frame.

Here you can observe the object internal structure. In this case, the `iris` object is a data frame with 5 variables. 
The first four variables are numeric, and the `iris$Species` is a factor with 3 levels, corresponding to the 3 Iris species.
This is important to check, because some functions only work with specific types of data.

For example, `ggplot(data =)` will expect a `data.frame` as input.

Other basic functions useful for data exploration:

```r
# View the first few rows
head(iris)

# Get a summary of all variables
summary(iris)

# Check the dimensions (rows x columns)
dim(iris)  # 150 rows, 5 columns

# See the three species
levels(iris$Species)
```

### Step 2: Compute Summary Statistics


```r
# Average petal length across all flowers
mean(iris$Petal.Length)

# Average petal length per species
tapply(iris$Petal.Length, iris$Species, mean)

# Correlation between petal length and petal width
cor(iris$Petal.Length, iris$Petal.Width)

# Histogram of petal length
hist(iris$Petal.Length)
```

You should find a strong positive correlation (~0.96) — flowers with longer petals also tend to have wider petals.

### Step 3: Visualize the Data

Install and load `ggplot2` if you haven't already:

```r
install.packages("ggplot2")
library(ggplot2)
```

Create a scatter plot of petal length vs. petal width, colored by species:

```r
ggplot(iris, aes(x = Petal.Length, y = Petal.Width, color = Species)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "Iris Petal Dimensions by Species",
    x = "Petal Length (cm)",
    y = "Petal Width (cm)"
  ) +
  theme_minimal()
```

The plot will show three well-separated clusters — *setosa* flowers are notably smaller, while *versicolor* and *virginica* overlap somewhat. This visual separation is why the iris dataset is so often used to demonstrate classification techniques.

### Step 4: Save Your Script

Save all of the above commands in a file called `iris_analysis.R` and run it with:

```bash
Rscript iris_analysis.R # if you use Linux/macOS
```

or open the script in Rstudio/Positron and click the `Source` button.

Congratulations — you've completed your first R data analysis!

---

## Troubleshooting

**R won't install on my machine.**
Make sure you're downloading the correct version for your operating system from [CRAN](https://cran.r-project.org/). On macOS, you may need to allow the installer in *System Settings → Privacy & Security*.

**`install.packages()` fails with a permissions error.**
Try running R as an administrator (Windows) or use `sudo` on Linux. Alternatively, install to a personal library by running `install.packages("pkg", lib = "~/R/library")`.

**A package loads but I get a "function not found" error.**
You may have forgotten to call `library(packagename)` at the start of your script. Installation and loading are two separate steps.

**My plot doesn't appear.**
In RStudio, plots appear in the Output pane. If using the terminal, call `dev.off()` after saving a plot to a file with `png()` or `pdf()`.

**I get a warning about package versions.**
Warnings (as opposed to errors) generally don't stop your code from running. Check whether your output looks correct — if it does, you can usually proceed. Update packages with `update.packages()`.

### 🚫 Common Beginner Mistakes

#### Forgetting Quotes

Wrong:

``` r
install.packages(ggplot2)
```

Correct:

``` r
install.packages("ggplot2")
```

#### Case Sensitivity

``` r
Mean(height)  # ❌
mean(height)  # ✅
```

R is case-sensitive.

---

## Getting Help

R has extensive built-in documentation. To look up any function, use `?` followed by the function name:

```r
?mean       # documentation for the mean() function
?ggplot     # documentation for ggplot()
help.start() # opens a browser-based help index
```

To search across all installed packages:

```r
??regression   # finds all help pages mentioning "regression"
```

Beyond the built-in docs, the R community is very active. Useful resources include:

- [Stack Overflow (R tag)](https://stackoverflow.com/questions/tagged/r) — for specific coding questions
- [R for Data Science (free book)](https://r4ds.hadley.nz/) — a thorough, beginner-friendly introduction
- [CRAN Task Views](https://cran.r-project.org/web/views/) — curated package lists by topic
- [Posit Community Forum](https://community.rstudio.com/) — for RStudio and tidyverse questions

Got a question? The most popular LLMs can handle everything from quick lookups to complex problems. 

Be aware of allucinations, though. Always double-check your code and test to see if the results are accurate.

- **Claude** (Anthropic) — [claude.ai](https://claude.ai)
- **ChatGPT** (OpenAI) — [chatgpt.com](https://chatgpt.com)
- **Gemini** (Google) — [gemini.google.com](https://gemini.google.com)
- **Copilot** (Microsoft) — [copilot.microsoft.com](https://copilot.microsoft.com)
- **Llama** (Meta, open-source) — [llama.meta.com](https://llama.meta.com)
- **Mistral** (Mistral AI, open-source) — [mistral.ai](https://mistral.ai)