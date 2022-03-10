## Step 1: Set working directory 
  # Put the baboon_data.csv file in your Desktop
  # In bottom right panel of RStudio, click Files > Home > Desktop > More > Set As Working Directory

## Load packages to analyze and plot data
install.packages("ggplot2")
install.packages("dplyr")
library(ggplot2)
library(dplyr)

## Load data
bab <- read.csv("baboon_data.csv")

## Preview data to ensure loaded successfully
View(bab)
str(bab) # structure of the data

## This dataset contains 1 row per fecal sample collected on an adult female.
## The column "female" identifies a female's unique ID (name).
## The column "cycle_day" identifies how many days before deturgescence the fecal sample was collected.
## The column "estrogen" identifies the female's estrogen concentration (based on the fecal sample). Units are nanograms estrogen per gram dried feces.
## The column "swelling_size" identifies how large a female's sexual swelling was on the day the fecal sample was collected.
## The column "alpha_consort" indicates whether or not (1=yes, 0=no) a female consorted (and likely mated) with an alpha male on the day the fecal sample was collected.
## The column "nonalpha_consort" indicates wehther or not (1=yes, 0=no) a female consorted (and likely mated) with a NON-alpha male on the day the fecal sample was collected.

## Modify a data type from "numeric" to "factor":
class(bab$female)                   # We want R to recognize unique females as a categorical variable instead of continuous
bab$female <- as.factor(bab$female) # This will adjust the class from integer (continuous) to factor (categories)
class(bab$female)                   # Ta-da!

## SAMPLE SIZES
# The line of code below will calculate the total number of rows in the dataset:
nrow(bab)
# The line of code below will calculate the number of unique values in the column "X".
# Replace "X" with the name of the correct column to calculate the number of unique females in the dataset.
# Make sure you don't put quotes around the name of the column!
length(unique(bab$female))

## DATA EXPLORATION
# We want to condense our dataset so that each row is one cycle day.
# For each cycle day, we want to calculate:
# 1. the total number of fecal samples for that cycle day,
# 2. mean and standard error of female sexual swelling size,
# 3. mean and standard error of female estrogen concentration,
# 4. the proportion of cases where the female consorted with an alpha male on that cycle day,
# and 5. the proportion of cases where the female consorted with a non-alpha male on that cycle day.
# We will use two functions to accomplish this: "group_by" and "summarize".
# For example, the code below will GROUP BY female, and SUMMARIZE her mean estrogen level across fecal samples.
bab %>%
  group_by(female) %>%
  summarize(estrogen_mean = mean(estrogen))
# Modify the code below to GROUP BY cycle day instead of female:
bab_new <-
  bab %>%
  group_by(cycle_day) %>% # modify this line so it GROUPS BY cycle_day, then run the code
  summarize(n = n(),                                                   # count the number of samples used to calculate each row
            swelling_size_mean = mean(swelling_size),                  # mean swelling size
            swelling_size_se = sd(swelling_size) / sqrt(n),            # standard error of swelling size
            estrogen_mean = mean(estrogen),                            # mean estrogen concentration
            estrogen_se = sd(estrogen) / sqrt(n),                      # standard error of estrogen concentration
            alpha_consorts = sum(alpha_consort),                       # count number of times an alpha male consorted a female on that cycle day
            alpha_consortship_probability = alpha_consorts/n,          # calculate the probability of alpha male consort
            nonalpha_consorts = sum(nonalpha_consort),                 # count number of times a non-alpha male consorted a female on that cycle day
            nonalpha_consortship_probability = nonalpha_consorts/n)    # calculate the probability of non-alpha male consort

# If you did this correctly, if you run the following line of code you should see a dataframe with "cycle_day" as the first column:
head(bab_new)
                      
## DATA VISUALIZATION: PART ONE
# Now that we've summarized our data, we can visualize how female sexual skin swelling sizes change over the course of her reproductive cycle.
# Run the code below:
ggplot(data = bab_new,
       aes(x = cycle_day,                                   # this indicates what the X-axis is
           y = swelling_size_mean,                          # this indicates what the Y-axis is
           ymin = swelling_size_mean - swelling_size_se,    # lower limit of the standard error bar
           ymax = swelling_size_mean + swelling_size_se)) + # upper limit of the standard error bar
  geom_point() +                                            # add points to the graph
  geom_errorbar() +                                         # add error bars to the graph
  geom_line() +                                             # connect points on the graph
  xlab("cycle day") +                                       # here you can change the x-axis label
  ylab("sex skin swelling size") +                        # here you can change the y-axis label
  scale_x_continuous(limits = c(-22, 12), breaks = seq(-20, 10, 5)) +
  theme_minimal()
# Save the plot as a PDF to your Desktop:
ggsave("swelling_size_plot.pdf", plot = last_plot(), device = "pdf", width = 5.5, height = 3.5, units = "in")
# Now you can open and view the plot on your Desktop.

# Now, let's modify the above code to make a new plot.
# Instead of plotting MEAN SWELLING SIZE on the y-axis, now we want to plot MEAN ESTROGEN on the y-axis.
# Try modifying the code above so that the y-axis is "estrogen_mean" instead of "swelling_size_mean".
# Don't forget to change the error bars as well on lines 75-76!
# You should also change the name of the y-axis label on line 81.
ggplot(data = bab_new,
       aes(x = cycle_day,                                   # this indicates what the X-axis is
           y = estrogen_mean,                          # this indicates what the Y-axis is
           ymin = estrogen_mean - estrogen_se,    # lower limit of the standard error bar
           ymax = estrogen_mean + estrogen_se)) + # upper limit of the standard error bar
  geom_point() +                                            # add points to the graph
  geom_errorbar() +                                         # add error bars to the graph
  geom_line() +                                             # connect points on the graph
  xlab("cycle day") +                                       # here you can change the x-axis label
  ylab("estrogen level") +                        # here you can change the y-axis label
  scale_x_continuous(limits = c(-22, 12), breaks = seq(-20, 10, 5)) +
  theme_minimal()
# Once your plot looks right, you can save it to your Desktop:
ggsave("estrogen_plot.pdf", plot = last_plot(), device = "pdf", width = 5.5, height = 3.5, units = "in")

## DATA VISUALIZATION: PART TWO
# We can also visualize the mating success of alpha males on different female cycle days.
# Run the code below:
ggplot(data = bab_new,
       aes(x = cycle_day,
           y = alpha_consortship_probability)) +
  geom_bar(stat = "identity", fill = "darkgray") +
  xlab("cycle day") +
  ylab("probability of consorting with alpha male") +
  scale_x_continuous(breaks = seq(-20,10,5)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_minimal()
# Save the plot as a PDF to your Desktop:
ggsave("alpha_consorts_plot.pdf", plot = last_plot(), device = "pdf", width = 5.5, height = 3.5, units = "in")

# Now, let's modify the above code to make a new plot.
# Instead of plotting the probability of a consorting with an ALPHA MALE on the y-axis,
#     we want to plot the probability of consorting with a NON-ALPHA MALE on the y-axis.
# Try modifying the code above to accomplish this.
# Don't forget to change the y-axis label as well on line 104.
ggplot(data = bab_new,
       aes(x = cycle_day,
           y = nonalpha_consortship_probability)) +
  geom_bar(stat = "identity", fill = "darkgray") +
  xlab("cycle day") +
  ylab("probability of consorting with nonalpha male") +
  scale_x_continuous(breaks = seq(-20,10,5)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_minimal()
# Once your plot looks right, you can save it to your Desktop:
ggsave("nonalpha_consorts_plot.pdf", plot = last_plot(), device = "pdf", width = 5.5, height = 3.5, units = "in")

# You're done!