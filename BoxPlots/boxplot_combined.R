# this will generate boxplots for body temperarture data of male red-sided garter snakes

#clear environment
rm(list=ls())+dev.off()

# Upload data
df <- read.table("2017_temps.tab", header = TRUE)
head(df)

summary(df)

# Run t-tests        
test1 <- t.test(x=df$BTvol,y=df$BTcourt)
test1

test2 <- t.test(x=df$Twarm,y=df$Tcool)
test2

test3 <- t.test(x=df$Thot,y=df$BTvol)
test3

test4 <- t.test(x=df$Thot,y=df$BTcourt)
test4

#run anova on behavioral max temps
#setup data
df_long <- data.frame(
  value = c(df$Thot, df$BTvol, df$BTcourt),
  group = factor(rep(c("Thot", "BTvol", "BTcourt"), each = nrow(df)))
)
#run anova
anova_result <- aov(value ~ group, data = df_long)
summary(anova_result)
# test what is different
TukeyHSD(anova_result)



# create Boxplot with y-axis limit to 50
boxplot(df, 
        xlab = "Thermal Condition",
        ylab = expression("Body Temperature ( "*degree*C*") "), 
        names = c(expression(CT[min]),
                  expression(T[cool]),
                  expression(T[warm]),
                  expression(T[hot]),
                  expression(BT[vol]), 
                  expression(BT[court]), 
                  expression(CT[max])),
        col = c("cornflowerblue", "#999933", "#999933", 
                "darkorange3", "orangered", "orangered", "red"),
        ylim = c(0, 46))



# Add significance bars and p-values manually (adjust positions as needed)


# test1 (bracket between boxes 5 and 6)
y_bracket1 <- 40  # height of the horizontal part of the bracket
tick_height1 <- 1  # height of the vertical ticks

# vertical ticks
segments(x0 = 5, x1 = 5, y0 = y_bracket1, y1 = y_bracket1 - tick_height1)
segments(x0 = 6, x1 = 6, y0 = y_bracket1, y1 = y_bracket1 - tick_height1)

# horizontal bar
segments(x0 = 5, x1 = 6, y0 = y_bracket1, y1 = y_bracket1)

# annotation
text(x = 5.5, y = y_bracket1 + 1, labels = "ns")


# test2 (bracket between boxes 2 and 3)
y_bracket2 <- 35
tick_height2 <- 1

segments(x0 = 2, x1 = 2, y0 = y_bracket2, y1 = y_bracket2 - tick_height2)
segments(x0 = 3, x1 = 3, y0 = y_bracket2, y1 = y_bracket2 - tick_height2)
segments(x0 = 2, x1 = 3, y0 = y_bracket2, y1 = y_bracket2)

text(x = 2.5, y = y_bracket2 + 1, labels = "****")

#test3

y_bracket3 <- 40
tick_height3 <- 1

segments(x0 = 5, x1 = 5, y0 = y_bracket3, y1 = y_bracket3 - tick_height3)
segments(x0 = 4, x1 = 4, y0 = y_bracket3, y1 = y_bracket3 - tick_height3)
segments(x0 = 4, x1 = 5, y0 = y_bracket3, y1 = y_bracket3)

text(x = 4.5, y = y_bracket3 + 1, labels = "*")

#test4

y_bracket4 <- 45
tick_height4 <- 2

segments(x0 = 6, x1 = 6, y0 = y_bracket4, y1 = y_bracket4 - tick_height4)
segments(x0 = 4, x1 = 4, y0 = y_bracket4, y1 = y_bracket4 - tick_height4)
segments(x0 = 4, x1 = 6, y0 = y_bracket4, y1 = y_bracket4)

text(x = 5, y = y_bracket4 + 1, labels = "***")

