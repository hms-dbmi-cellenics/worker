{{- define "imagePullSecret" }}
{{- printf "{\"auths\": {\"%s\": {\"auth\": \"%s\"}}}" .Values.imageCredentials.registry (printf "%s:%s" .Values.gitlab_deploy_user .Values.imageCredentials.gitlab_deploy_password | b64enc) | b64enc }}
{{- end }}