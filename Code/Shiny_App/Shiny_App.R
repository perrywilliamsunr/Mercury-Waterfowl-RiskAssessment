library(shiny)

# Load data
flank.hg <- c(150.9775, 1150.9775, 2150.9775, 3150.9775, 4150.9775, 5150.9775, 
              6150.9775, 7150.9775, 8150.9775, 9150.9775, 10150.9775, 11150.9775, 
              12150.9775, 13150.9775, 14150.9775, 15150.9775, 16150.9775, 17150.9775, 
              18150.9775, 19150.9775, 20150.9775, 21150.9775, 22150.9775, 23150.9775, 
              24150.9775, 25150.9775, 26150.9775, 27150.9775, 28150.9775, 29150.9775, 
              30150.9775, 31150.9775, 32150.9775, 33150.9775, 34150.9775, 35150.9775, 
              36150.9775, 37150.9775, 38150.9775, 39150.9775, 40150.9775, 41150.9775, 
              42150.9775, 43150.9775, 44150.9775, 45150.9775, 46150.9775, 47150.9775, 
              48150.9775)

wood.duck.probability <- c(0.125263004395559, 0.15694852405058, 0.205487157052684, 0.256546465845885, 
                           0.318167614524092, 0.391642189029852, 0.460554548674041, 0.550653084182239, 
                           0.622836072745453, 0.694664916776035, 0.762848155323626, 0.820365393829552, 
                           0.870258108868196, 0.901256171489282, 0.929421077849301, 0.944878445097181, 
                           0.960939941253672, 0.971126804574714, 0.9801470741412, 0.984625960877445, 
                           0.987688269483157, 0.99000062496094, 0.993021269503989, 0.993958710913901, 
                           0.995250296856447, 0.99591692185905, 0.997146011707602, 0.996854363268962, 
                           0.998000124992188, 0.998083453117513, 0.998604253900798, 0.999000062496094, 
                           0.999208382809408, 0.998812574214112, 0.998562589838135, 0.999187550778076, 
                           0.999062558590088, 0.998958398433431, 0.999291710934733, 0.999125054684082, 
                           0.999458367185384, 0.999354207028727, 0.999020894527425, 0.999458367185384, 
                           0.999645855467367, 0.999645855467367, 0.999458367185384, 0.999625023436035, 
                           0.999583359373373)

mallard.probability <- c(0.630398100118743, 0.685602983146887, 0.740912026331688, 0.797991792179655, 
                         0.852655042393184, 0.894048288648626, 0.935545695060725, 0.95925254671583, 
                         0.973751640522467, 0.985063433535404, 0.988438222611087, 0.990521425744224, 
                         0.990750578088869, 0.993000437472658, 0.991167218715497, 0.99045892965023, 
                         0.989521488240318, 0.988959023394371, 0.987417453075849, 0.985250921817386, 
                         0.984334312438806, 0.982167781180343, 0.978876320229986, 0.978543007728684, 
                         0.97681394912818, 0.974709913963711, 0.972439222548591, 0.97000187488282, 
                         0.967606191279712, 0.965877132679207, 0.96318980063746, 0.961898214694915, 
                         0.960752452971689, 0.959752515467783, 0.957252671708018, 0.957023519363373, 
                         0.95439868341562, 0.952065495906506, 0.951357206841239, 0.950836406057955, 
                         0.948857363081474, 0.947190800574964, 0.946836656042331, 0.946107534945733, 
                         0.943774347436619, 0.942753577901381, 0.943482698997979, 0.942107784930108, 
                         0.940024581796971)

# Define function to calculate probability
calculate_probability <- function(species, mercury) {
  if (species == "Wood ducks") {
    probability <- approxfun(flank.hg, wood.duck.probability, rule = 2)(mercury)
  } else {
    probability <- approxfun(flank.hg, mallard.probability, rule = 2)(mercury)
  }
  return(probability)
}

# Define UI for application
ui <- fluidPage(
  
  # App title
  titlePanel("Breast Tissue Mercury Probability Calculator"),
  
  # Sidebar layout
  sidebarLayout(
    sidebarPanel(
      
      # Input for species selection
      selectInput("species", "Select Species:",
                  choices = c("Wood ducks", "Mallards"),
                  selected = "Wood ducks"),
      
      # Input for flank feather mercury concentration
      numericInput("mercury", "Enter Flank Feather Mercury Concentration (ng/g):", value = 10000, min = 0, max = 50000)
    ),
    
    # Main panel for displaying output
    mainPanel(
      
      # Output: Display probability
      textOutput("probability"),
      
      # Output: Plot
      plotOutput("mercuryPlot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Calculate probability
  output$probability <- renderText({
    probability <- calculate_probability(input$species, input$mercury)
    paste("Probability of exceeding EPA threshold:", round(probability * 100, 2), "%")
  })
  
  # Generate plot
  output$mercuryPlot <- renderPlot({
    if (input$species == "Wood ducks") {
      species_prob <- wood.duck.probability
    } else {
      species_prob <- mallard.probability
    }
    interpolated_prob <- approxfun(flank.hg, species_prob, rule = 2)
    y_values <- interpolated_prob(flank.hg)  # Compute y-values using the interpolation function
    plot(flank.hg, y_values, type = "l", xlab = "Flank Feather Mercury Concentration (ng/g)", 
         ylab = "Probability of Exceeding Threshold", main = "Relationship between Mercury Concentration and Probability")
    abline(h = 0.3, col = "red", lty = 2)  # Add EPA threshold line
    points(input$mercury, calculate_probability(input$species, input$mercury), col = "blue", pch = 16)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
