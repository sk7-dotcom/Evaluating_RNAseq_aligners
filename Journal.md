2021/05/03

Now that I have finished the majority of the required part of this project I want to tabulate thoughts, concerns and things I still have no clue about. 
Firstly, this exercise was hard to do because it felt like there was never enough information to try and replicate the outcomes accuratly. 

I was not able to replicate the general upregulation of genes on day 7 partly becausse the HISAT2 and SALMON + KALLISTO never seemed to agree with eachother. 
Absolutly hated bringing together data with tximport and a regular coutn matrix. Even after making the reference level clear the data was flipped in the end. There has to be a trick, this is something I need to learn. 

I was not able to access day 14 data within contrast. All the pieces of this were just not available if I included an intercept. I should try to compare the data from a set with no intercept. Is there really a difference? What does it come down to? 

Don't really know why the HISAT2 system leads to greater sensitivity for low count gene calls. Need to look this up. 

What on earth have they done to the select function in dplyr. I litrally don't understand how it works anymore. Need to readjust code from previous work. 

What is the online vs offline phase of Salmon mean. Need to read more about this too. 

Reading about the HISAT2 mecahnism was the most fun I had in this project, lets learn more about Burrows Wheeler... I barely scratched the surface there. 

More thoughts to come, lets stay postive! 

2021/05/04 

Spent some more time thinking about next steps. Since there was a lot more I wanted to do in terms of finding out what pathways are actvated, why not drop HISAT2 for now and proceed with Salmon and Kallisto data. Start by fixing the treatment vs control problem and then ficuge out what the next steps are for pathway analysis. Not somehting I have done much reading on, will need a lot of help here from Omar. 

Also want to understand the difference between the Wald test and the LTR. Does that solve my issue with the correct fold changes? I have the code in the second R-script, maybe thats a good point of comparison. Since I got the day 14 data generated there, maybe thats a good way to compare day 7 to day 14 readily and see if I can go through the whole way. Perhaps I would have been better off taking this approach from the start. 
