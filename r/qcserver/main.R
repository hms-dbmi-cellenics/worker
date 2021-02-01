
source("qcserver/wrapper.R")

# Start remoter::server on its default port
# in debug mode for now (which might be desirable even in prod).
# The server will shut down when a client runs shutdown() in a session.
remoter::server(port = 55555, showmsg = TRUE)

