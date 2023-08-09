{{/* Generate a template that can be used for both assigned and unassigned experiments */}}
{{- define "worker.pod-template" -}}
    metadata:
      labels:
        type: 'worker'
        sandboxId: "{{ .Values.sandboxId }}"
    spec:
      containers:
      - name: "{{ .Release.Name }}-r"
        image: "{{ .Values.r.image.registry }}/{{ .Values.r.image.repository }}:{{ .Values.r.image.tag }}"
        volumeMounts:
        - name: 'data'
          mountPath: '/data'
        - name: watch-script
          mountPath: /var/lib/watchfile
          readOnly: true
        - name: shutdown-file
          mountPath: /var/lib/shutdown-file
        - name: podinfo
          mountPath: /etc/podinfo
        ports:
        - containerPort: 4000
        resources:
          requests:
            memory: "{{ .Values.r.memoryRequest }}"
      - name: "{{ .Release.Name }}"
        image: "{{ .Values.python.image.registry }}/{{ .Values.python.image.repository }}:{{ .Values.python.image.tag }}"
        env:
        - name: AWS_ACCOUNT_ID
          value: "{{ .Values.myAccount.accountId }}"
        - name: AWS_XRAY_DAEMON_ADDRESS
          value: xray-service.default:2000
        - name: 'K8S_ENV'
          value: {{ .Values.kubernetes.env | quote }}
        - name: 'IGNORE_TIMEOUT'
          valueFrom:
            configMapKeyRef:
              name: instance-config
              key: ignoreTimeout
        volumeMounts:
        - name: 'data'
          mountPath: '/data'
        - name: watch-script
          mountPath: /var/lib/watchfile
          readOnly: true
        - name: shutdown-file
          mountPath: /var/lib/shutdown-file
        - name: podinfo
          mountPath: /etc/podinfo
        resources:
          requests:
            memory: "1Gi"
{{- if eq .Values.clusterEnv "production" }}
      - name: datadog-agent
        image: datadog/agent
        env:
        - name: DD_API_KEY
          value: "{{ .Values.myAccount.datadogApiKey }}"
        - name: DD_SITE
          value: "datadoghq.eu"
        - name: DD_EKS_FARGATE
          value: "true"
        - name: DD_CLUSTER_NAME
          value: "biomage-{{ .Values.kubernetes.env }}"
        - name: DD_TAGS
          value: "{{ .Values.datadogTags }}"
        - name: DD_KUBERNETES_POD_LABELS_AS_TAGS
          value: '{"*": "%%label%%"}'
        # Disable log collection by DD agent
        # because we push logs to Cloudwatch
        - name: DD_LOGS_ENABLED
          value: "false"
        - name: DD_CONTAINER_EXCLUDE
          value: "name:.*"
        - name: DD_CONTAINER_INCLUDE_METRICS
          value: "name:worker name:worker-r"
        - name: DD_KUBERNETES_KUBELET_NODENAME
          valueFrom:
            fieldRef:
              apiVersion: v1
              fieldPath: spec.nodeName
{{- end }}
      volumes:
      - name: 'data'
      - name: watch-script
        configMap:
          name: "watch-script"
          items:
            - key: watcher.sh
              path: watcher.sh
            - key: entrypoint.sh
              path: entrypoint.sh
      - name: shutdown-file
      - name: podinfo
        downwardAPI:
          items:
            - path: "labels"
              fieldRef:
                fieldPath: metadata.labels
      restartPolicy: Always
      serviceAccountName: 'deployment-runner'
{{- end -}}
