#Number of people of a certain age separated by sex. Males and females are positioned one at the side of each other.
data_rancagua <- read.csv("DataRancagua.csv")
data_all[, 'Date'] <- as.Date(data_all[, 'Date'])
mygraph = ggplot(data_all, aes(Date,Age))
mygraph + geom_histogram(aes(fill=Sex, color=Sex), position = 'dodge') + theme_minimal() + scale_x_continuous(breaks=seq(0,80,by=4)) + labs(x="Age (Years)", y="Frequency")
ggsave("FigSup2A.png")

#Number of both females and males tested in each month - line graph.
data_rancagua <- read.csv("TimelineRancaguaBySex_2B.csv")
ggplot(data_rancagua, aes(Month,Population, group=Sex, color=Sex)) + geom_line() + geom_point() + theme_minimal() + theme(text = element_text(size=15)) +theme(axis.text.x = element_text(angle=90, hjust=1))+ labs(x="Month", y="Age (Years)")
ggsave("FigSup2B.png")

#Number of both females and males tested in each month - bar graph.
ggplot(data_rancagua, aes(Population, Month, fill=Sex)) + geom_bar(stat ="identity", position="dodge") + coord_flip() + theme_minimal() + theme(axis.text.x = element_text(angle=90, hjust=1))+ labs(x="Population", y="Month")
ggsave("FigSup2A_2.png")

#Number of both females and males by age range (decade) - bar graph with males and females positioned one at the side of each other.
data_rancagua <- read.csv("TimelineRancaguaByAge.csv")
ggplot(data_rancagua, aes(Age.Range, fill=Sex)) + geom_bar(position="dodge") + theme_minimal() + theme(axis.text.x = element_text(angle=90, hjust=1))+ labs(x="Age Range (Years)", y="Population") 
ggsave("SupFig2E.png")

#Number of tested people for age range (decade) - bar graph with age ranges pilled up.
ggplot(data_rancagua, aes(Month)) + geom_bar(aes(fill= AgeRange)) + theme_minimal() + theme(axis.text.x = element_text(angle=90, hjust=1))+ labs(x="Month", y="Population")
ggsave("SupFig2C.png")

#Number of tested people of each age range (decade) per month - line graph.
ggplot(data_rancagua, aes(Month,Population, group=AgeRange, color=AgeRange)) + geom_line() + geom_point() + theme_minimal() + theme(text = element_text(size=15)) +theme(axis.text.x = element_text(angle=90, hjust=1))+ labs(x="Month", y="Population") 
ggsave("SupFig2D.png")
