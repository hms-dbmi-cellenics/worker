library(RestRserve)

app <- Application$new()

app$add_get(
    path = "/hello",
    FUN = function(request, response) {
        response$set_body("up")
    }
)

backend <- BackendRserve$new()
backend$start(app, http_port = 4000)