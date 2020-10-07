library(RestRserve)

app <- Application$new()

app$add_get(
    path = "/hello",
    FUN = function(request, response) {
        response$set_body("Hello from RestRserve")
    }
)

backend <- BackendRserve$new()
backend$start(app, http_port = 8080)